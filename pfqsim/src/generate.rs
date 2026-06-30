//! n2bio/pfqsim/src/generate.rs
//! 

use std::io::{self, Error, ErrorKind};
use std::thread;
use crossbeam_channel::bounded;
use rayon::prelude::*;
use rand::rngs::SmallRng;
use std::time::Instant;

use n2core::fasta::FastaReader;
use n2core::readers::ReaderType;
use n2core::writers::WriterType;
use n2core::fastq::{FastqRecord, PairedFastqRecord, Read1, Read2, PairedFastqWriter };
use n2core::sequence::DnaSequence; 

use crate::cli::GenerateArgs;
use crate::simstats::LibraryModel;
use crate::inserts::InsertSize;
use crate::qualities::QualityScores;
use crate::genome::ReferenceGenome;
use crate::mutate::{Mutator, MutationStats};

// ============================================================================
// Main Runner
// ============================================================================

pub(crate) fn run(args: GenerateArgs) -> io::Result<()> {
    let start_time: Instant = Instant::now();
    println!("Generating {} reads from {:?}", args.num_reads, args.fasta);

    let read_length: usize = args.length; 
    let deletion_buffer: usize = 20; 
    
    // Assume a max insert size of ~1000, plus the buffer
    let min_required_length: usize = 1000 + deletion_buffer;
   
    // 1. Initialize the FastaReader and load the weighted ReferenceGenome
    let ref_reader: FastaReader<ReaderType> = FastaReader::open(args.fasta.to_str().unwrap())?;
    let reference: ReferenceGenome = ReferenceGenome::load(ref_reader, min_required_length, args.circular)?;

    // 2. Load model, initialize mutator and samplers
    let model_path: &str = args.model.to_str().ok_or_else(|| {
        Error::new(ErrorKind::InvalidInput, "Model path is not valid UTF-8")
    })?;
    
    let model: LibraryModel = LibraryModel::from_file(model_path)
        .map_err(|e| Error::new(ErrorKind::Other, e.to_string()))?;
    
    // 3. Instantiate components
    let inserts: InsertSize = InsertSize::new(&model.insert_size)
        .map_err(|e| Error::new(ErrorKind::Other, e))?;
    let qualities: QualityScores = QualityScores::new(&model.quality)
        .map_err(|e| Error::new(ErrorKind::Other, e))?;
    
    let mutator = Mutator::new(args.sub_rate, args.indel_rate);

    // 4. Setup global thread pool for Rayon based on user args
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap_or_else(|_| ());

    // 5. Create a SINGLE bounded channel for batches of PairedReads
    let (tx, rx) = bounded::<Vec<PairedFastqRecord>>(50);

    // 6. Spawn the dedicated Writer Thread
    let prefix: String = args.prefix.clone();
    let gz_threads: usize = if args.threads > 2 { 2 } else { 1 };

    let writer_handle: thread::JoinHandle<Result<usize, Error>> = thread::spawn(move || -> io::Result<usize> {
        let r1_writer: WriterType = WriterType::to_multithreaded_gz(&format!("{}.r1.fq.gz", prefix), gz_threads)?;
        let r2_writer: WriterType = WriterType::to_multithreaded_gz(&format!("{}.r2.fq.gz", prefix), gz_threads)?;
        let mut fastq_writer: PairedFastqWriter<WriterType, WriterType> = PairedFastqWriter::new(r1_writer, r2_writer);
        
        let mut total_pairs_written: usize = 0;

        for read_batch in rx {
            for pair in read_batch {
                fastq_writer.write_pair(&pair)?; 
                total_pairs_written += 1;
            }
        }
        
        Ok(total_pairs_written)
    });

    // 7. Worker Pool (Rayon) chunking logic
    let batch_size: usize = 10_000;
    let num_batches: usize = (args.num_reads as f64 / batch_size as f64).ceil() as usize;

    // 8. Set min insert size
    let min_insert_size: usize = read_length + deletion_buffer;

    // 9. Iterate over batches
    (0..num_batches).into_par_iter().for_each_with(tx, |sender, batch_idx| {
        
        let mut rng: SmallRng = rand::make_rng();
        let mut batch: Vec<PairedFastqRecord> = Vec::with_capacity(batch_size);
        
        let reads_this_batch: usize = std::cmp::min(batch_size, args.num_reads - (batch_idx * batch_size));
        let start_read_id: usize = batch_idx * batch_size;

        for local_i in 0..reads_this_batch {
            let global_read_id: usize = start_read_id + local_i;

            // A. Sample an insert size
            let mut insert_size: usize = inserts.sample(&mut rng);
            while insert_size < min_insert_size {
                insert_size = inserts.sample(&mut rng);
            }

            // B. Grab a slice of the genome equal to the INSERT SIZE
            let (accession, raw_insert_slice) = reference.sample_slice(&mut rng, insert_size, deletion_buffer);

            // C. Mutate R1 from the beginning of the insert slice
            let r1_stats: MutationStats = mutator.mutate(raw_insert_slice, read_length, &mut rng);

            // D. Mutate R2 from the end of the insert slice
            let r2_start: usize = raw_insert_slice.len().saturating_sub(read_length + deletion_buffer);
            let mut r2_stats: MutationStats = mutator.mutate(&raw_insert_slice[r2_start..], read_length, &mut rng);
            
            // E. Reverse complement Read 2
            r2_stats.sequence = r2_stats.sequence.reverse_complement();

            // F. Format headers
            let r1_base_id: String = format!("{}:{}:{} 1:N:{}:{}:{}",
                args.prefix, accession, global_read_id, r1_stats.subs, r1_stats.insertions, r1_stats.deletions, 
            );

            let r2_base_id: String = format!("{}:{}:{} 2:N:{}:{}:{}",
                args.prefix, accession, global_read_id, r2_stats.subs, r2_stats.insertions, r2_stats.deletions, 
            );

            // G. Generate qualities
            let r1_qual: Vec<u8> = qualities.generate(&mut rng, read_length, 1);
            let r2_qual: Vec<u8> = qualities.generate(&mut rng, read_length, 2);

            // H. Convert Vec<u8> buffers to Strings
            let (r1_seq_str, r1_qual_str, r2_seq_str, r2_qual_str) = unsafe {
                (
                    String::from_utf8_unchecked(r1_stats.sequence),
                    String::from_utf8_unchecked(r1_qual),
                    String::from_utf8_unchecked(r2_stats.sequence),
                    String::from_utf8_unchecked(r2_qual),
                )
            };

            // I. Package and push to batch
            let r1_record: FastqRecord<Read1> = FastqRecord::<Read1>::new(r1_base_id, r1_seq_str, r1_qual_str);
            let r2_record: FastqRecord<Read2> = FastqRecord::<Read2>::new(r2_base_id, r2_seq_str, r2_qual_str);

            batch.push(PairedFastqRecord {
                r1: r1_record,
                r2: r2_record,
            });
        }

        sender.send(batch).expect("Failed to send batch to writer");
    });
    
    // 10. Run summary
    let pairs_written: usize = writer_handle.join().expect("Writer thread panicked")?;
    let duration: std::time::Duration = start_time.elapsed();
    let summary: serde_json::Value = serde_json::json!({
        "source_fasta": args.fasta,
        "total_pairs_generated": pairs_written,
        "runtime_seconds": duration.as_secs_f64()
    });

    println!("{}", serde_json::to_string_pretty(&summary).unwrap());

    Ok(())
}