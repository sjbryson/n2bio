//! n2bio/pfqbz2gz/src/main.rs

use std::io;
use std::sync::mpsc;
use std::thread;
use std::time::Instant;
use clap::Parser;
use n2core::fastq::{PairedFastqReader, PairedFastqWriter, PairedFastqRecord};
use n2core::writers::WriterType; 

#[derive(Parser, Debug)]
#[command(author, version, about = "Convert paired bz2 FASTQ to gz FASTQ using multi-threaded pipeline.")]
struct Args {
    /// Path to R1 bz2 file
    #[arg(short = '1', long = "r1")]
    r1: String,

    /// Path to R2 bz2 file
    #[arg(short = '2', long = "r2")]
    r2: String,

    /// Output prefix for the new gz files (e.g. 'sample1' becomes 'sample1_R1.fq.gz')
    #[arg(short = 'o', long = "output-prefix")]
    output_prefix: String,

    /// Total CPU threads to allocate across the pipeline
    #[arg(short = 't', long = "threads", default_value_t = 4)]
    threads: usize,
}

fn main() -> io::Result<()> {
    let start_time = Instant::now();
    let args: Args = Args::parse();
    let out_r1: String = format!("{}_R1.fq.gz", args.output_prefix);
    let out_r2: String = format!("{}_R2.fq.gz", args.output_prefix);
    let gz_threads: usize = if args.threads > 2 { 2 } else { 1 };
    // Create a bounded channel for the reader.
    let (tx, rx) = mpsc::sync_channel::<PairedFastqRecord>(2000);

    // Spawn the Compression / Writer Thread
    let writer_handle: thread::JoinHandle<Result<usize, io::Error>> = thread::spawn(move || -> io::Result<usize> {
        let r1_writer: WriterType = WriterType::to_multithreaded_gz(&out_r1, gz_threads)?;
        let r2_writer: WriterType = WriterType::to_multithreaded_gz(&out_r2, gz_threads)?;
        
        let mut fastq_writer: PairedFastqWriter<WriterType, WriterType> = PairedFastqWriter::new(r1_writer, r2_writer);
        let mut pairs_written: usize = 0;

        // Block and wait for decompressed records from the main thread
        for pair in rx {
            fastq_writer.write_pair(&pair)?;
            pairs_written += 1;
        }
        
        fastq_writer.flush()?;
        Ok(pairs_written)
    });

    // Decompression / Reader Loop
    let paired_reader: PairedFastqReader<n2core::readers::ReaderType, n2core::readers::ReaderType> = PairedFastqReader::from_bzs(&args.r1, &args.r2)?;
    for result in paired_reader {

        // Send to the writer thread. Break the loop if writer crashes.
        if tx.send(result?).is_err() {
            eprintln!("Error: Writer thread disconnected unexpectedly.");
            break;
        }
    }

    // Clean up and wait
    drop(tx);

    // Wait for writer to finish.
    let total_records: usize = writer_handle.join().expect("Writer thread panicked")?;

    let duration: std::time::Duration = start_time.elapsed();

    let summary = serde_json::json!({
        "input_R1"        : args.r1,
        "input_R2"        : args.r2,
        "output_prefix"   : args.output_prefix,
        "threads"         : args.threads,
        "total_pairs"     : total_records,
        "runtime_seconds" : duration.as_secs_f32(),
    });

    println!("{}", serde_json::to_string_pretty(&summary).unwrap());

    Ok(())
}