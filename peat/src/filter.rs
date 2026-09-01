//! n2bio/peat/src/filter.rs
//! 

use clap::Parser;
use crossbeam::channel::bounded;
use std::io::{self, BufRead, Write};
use std::thread;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

use n2bio::sam::{SamReader, SamStr, SamFields, SamFlags, SamTags, AlignmentStats};
use n2bio::fastq::{ShardedMateMap, PairedRead, PairedFastqWriter};
use n2bio::writers::WriterType;

use crate::cli::{ FilterArgs, ThresholdMode, ThresholdMetrics };
use crate::samfilters::{ lowpass_filter, highpass_filter, threshold_args };



pub(crate) fn run(args: FilterArgs) -> io::Result<()> {

    let start_time: Instant = Instant::now(); // Start the clock

    // Define filter function based on ThresholdMode (lowpass or highpass) and optional thresholds 
    let has_thresholds: bool = threshold_args(&args.thresholds);
    let sam_filter: fn(&SamStr<'_>, &ThresholdMetrics, bool) -> bool = match args.filter_mode {
        ThresholdMode::LowPass => lowpass_filter,
        ThresholdMode::HighPass => highpass_filter,
    };
    
    // Initialize shared state & channels.
    let mate_map: Arc<ShardedMateMap> = Arc::new(ShardedMateMap::new(args.shards));
    let total_primary_reads: Arc<AtomicU64> = Arc::new(AtomicU64::new(0));
    
    // Channels act as elastic buffers between pipeline stages.
    let (line_tx, line_rx) = bounded::<String>(10_000);
    let (pair_tx, pair_rx) = bounded(10_000);

    // Spawn FASTQ writer thread using multithreaded gzip.
    // Dedicate a couple of internal threads to the gzp compression engine.
    let fq_prefix: String = args.prefix.clone();
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
    let mut worker_handles: Vec<thread::JoinHandle<()>> = Vec::with_capacity(args.threads);
    for _ in 0..args.threads {
        let rx: crossbeam::channel::Receiver<String> = line_rx.clone();
        let p_tx: crossbeam::channel::Sender<n2core::fastq::PairedFastqRecord> = pair_tx.clone();
        let map: Arc<ShardedMateMap> = Arc::clone(&mate_map);
        let primary_counter: Arc<AtomicU64> = Arc::clone(&total_primary_reads);
        let worker_args: FilterArgs = args.clone();

        let handle: thread::JoinHandle<()> = thread::spawn(move || {
            for line in rx {               
                let sam: SamStr<'_> = SamStr::new(&line);
                if !sam.is_primary() { continue; }
                primary_counter.fetch_add(1, Ordering::Relaxed);
                
                if sam_filter(&sam, &worker_args.thresholds, has_thresholds) {

                    // If passes filter - parse strictly and push to the ShardedMateMap.
                    let paired_read: PairedRead = PairedRead::from_samstr(&sam);
                        
                    // If the pair resolves - forward it to the FASTQ writer thread.
                    if let Some(pair) = map.process(paired_read) {
                        let _ = p_tx.send(pair); 
                        // Update logic here for "--stdout interleaved" feature
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
        // Need to update logic for "--stdout sam" feature
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

    // Update to use write to file feature "--report {args.prefix}.filter_report.json"
    // Summary data.
    let duration: std::time::Duration = start_time.elapsed();
    let primary_reads: u64 = total_primary_reads.load(Ordering::Relaxed);
    let total_pairs: u64 = primary_reads / 2;
    
    let summary: serde_json::Value = serde_json::json!({
        "runtime_seconds": duration.as_secs_f64(),
        "total_pairs"    : total_pairs,
        "written_pairs"  : pairs_written,
        "orphaned_reads" : mate_map.orphan_count(),
        "AP"             : args.thresholds.align_prop.map_or(serde_json::Value::Null, |v| serde_json::json!(v)),
        "AI"             : args.thresholds.align_ident.map_or(serde_json::Value::Null, |v| serde_json::json!(v)),
        "AS"             : args.thresholds.align_score.map_or(serde_json::Value::Null, |v| serde_json::json!(v)),
        "AL"             : args.thresholds.align_length.map_or(serde_json::Value::Null, |v| serde_json::json!(v)),
        "BS"             : args.thresholds.base_score.map_or(serde_json::Value::Null, |v| serde_json::json!(v)),
        "MQ"             : args.thresholds.mapq.map_or(serde_json::Value::Null, |v| serde_json::json!(v)),
    });

    // Write JSON to file - {args.report}.json
    let out_json: String = format!("{}.json", args.report);
    let mut json_file: std::fs::File = std::fs::File::create(&out_json)?;
    json_file.write_all(serde_json::to_string_pretty(&summary).unwrap().as_bytes())?;

    Ok(())
}

