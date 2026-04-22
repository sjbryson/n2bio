//! n2bio/fastfilter/src/main.rs

use clap::Parser;
use crossbeam::channel::bounded;
use std::io::{self, BufRead};
use std::thread;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;
use n2core::sam::{SamReader, SamStr, SamFields, SamFlags, SamTags, AlignmentStats};
use n2core::fastq::{ShardedMateMap, PairedRead, PairedFastqWriter};
use n2core::writers::WriterType;


#[derive(Parser, Debug, Clone)]
#[command(author, version, about = "High-performance SAM stream to paired FASTQ filter", long_about = None)]
struct Args {

    /// Number of worker threads for parsing and pairing
    #[arg(short = 't', long, default_value_t = 4)]
    threads: usize,

    /// Number of shards for the ShardedMateMap (recommend 4-8x threads)
    #[arg(long, default_value_t = 64)]
    shards: usize,

    /// Prefix for output files (e.g. 'out' -> out.r1.fq.gz, out.r2.fq.gz)
    #[arg(short = 'p', long, required = true)]
    fq_prefix: String,

    /// If none of the following optional max stats are set, only unmapped pairs are written.
    /// Optional: Maximum Alignment Proportion - sam.calculate_alignment_proportion()
    #[arg(long)]
    max_ap: Option<f32>,
    
    /// Optional: Maximum Percent Identity - sam.calculate_alignment_accuracy()
    #[arg(long)]
    max_pi: Option<f32>,
    
    /// Optional: Max AS score for filtering mapped reads - sam.get_int_tag("AS")
    #[arg(long)]
    max_as: Option<i32>,

    /// Optional: Maximum Alignment Lenth - sam.calculate_alignment_length()
    #[arg(long)]
    max_al: Option<u32>,
    
    /// Optional: Max AS/AL score for filtering mapped reads - sam.calculate_as_al()
    #[arg(long)]
    max_sl: Option<f32>,

    /// Optional: Max MAPQ score for filtering mapped reads - sam.calculate_as_al()
    #[arg(long)]
    max_mq: Option<u32>,
}

/// Filter logic for whether an alignment passes - didn't align well in this use case.
fn sam_filter(sam: &SamStr, args: &Args) -> bool {
    // Keep unmapped reads
    if !sam.is_mapped() {
        return true;
    }
    // Evaluate optional filters
    if args.max_ap.is_some_and(|max: f32| sam.calculate_alignment_proportion().is_some_and(|val: f32| val <= max)) {
        return true;
    }
    if args.max_pi.is_some_and(|max: f32| sam.calculate_alignment_accuracy().is_some_and(|val: f32| val <= max)) {
        return true;
    }
    if args.max_as.is_some_and(|max: i32| sam.get_int_tag("AS").is_some_and(|val: i32| val <= max)) {
        return true;
    }
    if args.max_al.is_some_and(|max: u32| sam.calculate_alignment_length().is_some_and(|val: u32| val <= max)) {
        return true;
    }
    if args.max_sl.is_some_and(|max: f32| sam.calculate_as_al().is_some_and(|val: f32| val <= max)) {
        return true;
    }
    if args.max_mq.is_some_and(|max: u32| sam.mapq() <= max) {
        return true;
    }
    // If it is mapped but didn't pass any of the specified max thresholds
    false
}

fn main() -> io::Result<()> {

    let start_time: Instant = Instant::now(); // Start the clock!
    let args: Args = Args::parse();
    
    // Initialize shared state & channels.
    let mate_map: Arc<ShardedMateMap> = Arc::new(ShardedMateMap::new(args.shards));
    let total_primary_reads: Arc<AtomicU64> = Arc::new(AtomicU64::new(0));
    
    // Channels act as elastic buffers between pipeline stages.
    let (line_tx, line_rx) = bounded::<String>(10_000);
    let (pair_tx, pair_rx) = bounded(10_000);

    // Spawn FASTQ writer thread using multithreaded gzip.
    // Dedicate a couple of internal threads to the gzp compression engine.
    let fq_prefix: String = args.fq_prefix.clone();
    let gz_threads: usize = if args.threads > 2 { 2 } else { 1 }; 
    let fastq_writer_handle: thread::JoinHandle<Result<usize, io::Error>> = thread::spawn(move || -> io::Result<usize> {
        let r1_writer: WriterType = WriterType::to_multithreaded_gz(&format!("{}.r1.fq.gz", fq_prefix), gz_threads)?;
        let r2_writer: WriterType = WriterType::to_multithreaded_gz(&format!("{}.r2.fq.gz", fq_prefix), gz_threads)?;
        let mut fastq_writer: PairedFastqWriter<WriterType, WriterType> = PairedFastqWriter::new(r1_writer, r2_writer);
        let mut pairs_written: usize = 0;

        for pair in pair_rx {
            fastq_writer.write_pair(&pair)?;
            pairs_written += 1;
        }
        
        fastq_writer.flush()?;
        Ok(pairs_written)
    });

    // Spawn worker threads.
    let mut worker_handles = Vec::with_capacity(args.threads);
    for _ in 0..args.threads {
        let rx: crossbeam::channel::Receiver<String> = line_rx.clone();
        let p_tx: crossbeam::channel::Sender<n2core::fastq::PairedFastqRecord> = pair_tx.clone();
        let map: Arc<ShardedMateMap> = Arc::clone(&mate_map);
        let primary_counter: Arc<AtomicU64> = Arc::clone(&total_primary_reads);
        let worker_args: Args = args.clone();

        let handle: thread::JoinHandle<()> = thread::spawn(move || {
            for line in rx {               
                let sam: SamStr<'_> = SamStr::new(&line);
                if !sam.is_primary() { continue; }
                primary_counter.fetch_add(1, Ordering::Relaxed);
                
                if sam_filter(&sam, &worker_args) {

                    // If passes filter - parse strictly and push to the ShardedMateMap.
                    let paired_read: PairedRead = PairedRead::from_samstr(&sam);
                        
                    // If the pair resolves - forward it to the FASTQ writer thread.
                    if let Some(pair) = map.process(paired_read) {
                        let _ = p_tx.send(pair);
                    }
                }
            }
        });
        worker_handles.push(handle);
    }

    // Drop main thread's copies of the senders.
    // This allows channels to close when workers finish.
    drop(pair_tx);
   
    // Setup the Reader on the main thread.
    let mut sam_reader: SamReader = SamReader::from_stdin();

    // The streaming loop.
    let mut line_buffer: String = String::new();
    while let Ok(bytes) = sam_reader.reader.read_line(&mut line_buffer) {
        if bytes == 0 { break; } // Clean EOF.
        let clean_line: String = line_buffer.trim_end().to_string();
        // Skip SAM header lines
        if clean_line.starts_with('@') {
            line_buffer.clear();
            continue;
        }
        // Push to workers. If the queue is full, the main thread safely blocks.
        if line_tx.send(clean_line).is_err() {
            break; 
        }
        line_buffer.clear();
    }

    // Tear down.
    drop(line_tx); // Signals workers that the file is fully read.

    // Wait for workers to finish mapping.
    for handle in worker_handles {
        handle.join().unwrap();
    }

    // Wait for writers to finish flushing to disk.
    let pairs_written: usize = fastq_writer_handle.join().unwrap().unwrap();

    // Summary data.
    let duration: std::time::Duration = start_time.elapsed();
    let primary_reads: u64 = total_primary_reads.load(Ordering::Relaxed);
    let total_pairs: u64 = primary_reads / 2;
    
    let summary = serde_json::json!({
        "runtime_seconds": duration.as_secs_f64(),
        "total_pairs"    : total_pairs,
        "written_pairs"  : pairs_written,
        "orphaned_reads" : mate_map.orphan_count(),
        "max_ap"         : args.max_ap.map_or(serde_json::Value::Null, |v| serde_json::json!(v)),
        "max_pi"         : args.max_pi.map_or(serde_json::Value::Null, |v| serde_json::json!(v)),
        "max_as"         : args.max_as.map_or(serde_json::Value::Null, |v| serde_json::json!(v)),
        "max_al"         : args.max_al.map_or(serde_json::Value::Null, |v| serde_json::json!(v)),
        "max_sl"         : args.max_sl.map_or(serde_json::Value::Null, |v| serde_json::json!(v)),
        "max_mq"         : args.max_mq.map_or(serde_json::Value::Null, |v| serde_json::json!(v)),
    });

    println!("{}", serde_json::to_string_pretty(&summary).unwrap());

    Ok(())
}