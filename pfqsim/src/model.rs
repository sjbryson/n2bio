//! n2bio/pfqsim/src/model.rs
//! 

use std::io::{ self, Error, ErrorKind };
use std::time::Instant;

use n2core::bam::{ BamReader, BamRecord, BamStats, BamFlags };

use crate::cli::ModelArgs;
use crate::simstats::{LibraryModel, LibraryProfiler};

// ============================================================================
// Main Runner
// ============================================================================

pub(crate) fn run(args: ModelArgs) -> io::Result<()> {
    let start_time: Instant              = Instant::now();
    let read_length: usize               = args.length;
    let mut total_pairs: usize           = 0;
    let max_insert_size: i32             = args.max_ins as i32;
    let mut bam_reader: BamReader        = BamReader::open(args.bam.to_str().unwrap())?;
    let _header: n2core::bam::BamHeader  = bam_reader.read_header()?;
    let mut current_record: BamRecord    = BamRecord::default();
    let mut r1_record: Option<BamRecord> = None;
    let mut r2_record: Option<BamRecord> = None;
    let mut prev_qname: Vec<u8>          = Vec::new();
    let mut profiler: LibraryProfiler    = LibraryProfiler::new();

    let mut process_pair = |r1: &BamRecord, r2: &BamRecord| {
        // Accumulate base quality distributions
        profiler.add_qualities(&r1.qual, &r2.qual);
        // Filter and accumulate insert sizes
        if r1.mapq as usize >= args.mapq && r2.mapq as usize >= args.mapq {
            if r1.ref_id == r2.ref_id && r1.ref_id != -1 {
                let (fwd, rev) = if r1.pos <= r2.pos {
                    (r1, r2)
                } else {
                    (r2, r1)
                };
                let ref_span: i32 = rev.calculate_ref_span().unwrap_or(0) as i32;
                let insert_size: i32 = (rev.pos + ref_span) - fwd.pos;
                if insert_size > 0 && insert_size <= max_insert_size {
                    profiler.add_insert_size(insert_size);
                }
            }
        }
    };

    // Main name-sorted parsing loop
    while bam_reader.read_record(&mut current_record)? {
        if current_record.read_name != prev_qname {
            if !prev_qname.is_empty() {
                total_pairs += 1;
                if let (Some(r1), Some(r2)) = (&r1_record, &r2_record) {
                    process_pair(r1, r2);
                }
            }
            r1_record = None;
            r2_record = None;
            prev_qname = current_record.read_name.clone();
        }

        if current_record.is_read1() {
            r1_record = Some(current_record.clone());
        } else if current_record.is_read2() {
            r2_record = Some(current_record.clone());
        }
    }

    // Parse the last record pair in the file
    if !prev_qname.is_empty() {
        total_pairs += 1;
        if let (Some(r1), Some(r2)) = (&r1_record, &r2_record) {
            process_pair(r1, r2);
        }
    }

    // Check validity using the quality profiler
    if total_pairs > 0 && profiler.quality.r1_counts.is_empty() {
        return Err(Error::new(
            ErrorKind::InvalidData,
            "No valid pairs found. Is the BAM file name-sorted (`samtools sort -n`)?",
        ));
    }
    let insert_sizes_modeled: usize = profiler.insert_size.total_modeled;
    // Compile the final nested profile
    let final_model: LibraryModel = profiler.build(read_length);
    // Write model to file
    let out_path: &str = args.model.to_str().ok_or_else(|| {
        Error::new(ErrorKind::InvalidInput, "Output path is not valid UTF-8")
    })?;
    
    final_model.to_file(out_path).map_err(|e| {
        Error::new(ErrorKind::Other, format!("Failed to save JSON model: {}", e))
    })?;

    // 4. Output a summary to stdout
    let duration: std::time::Duration = start_time.elapsed();
    let summary: serde_json::Value = serde_json::json!({
        "total_pairs_scanned":  total_pairs,
        "insert_sizes_modeled": insert_sizes_modeled,
        "runtime_seconds":      duration.as_secs_f64()
    });

    println!("{}", serde_json::to_string_pretty(&summary).unwrap());

    Ok(())
}