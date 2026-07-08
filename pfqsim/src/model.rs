//! n2bio/pfqsim/src/model.rs
//! 

use std::io::{ self, Error, ErrorKind };
use std::time::Instant;
use std::fs::File;

use n2core::bam::{ BamReader, BamRecord, BamStats, BamFlags };

use crate::cli::ModelArgs;
use crate::modelstats::ModelStats;

// ============================================================================
// Main Runner
// ============================================================================

pub(crate) fn run(args: ModelArgs) -> io::Result<()> {
    let start_time: Instant              = Instant::now();
    let read_length: usize               = args.read_length;
    let mut total_pairs: usize           = 0;
    let mut insert_sizes_modeled: usize  = 0;
    let max_insert_size: i32             = args.max_ins as i32;
    let mut bam_reader: BamReader        = BamReader::open(&args.bam)?;
    let _header: n2core::bam::BamHeader  = bam_reader.read_header()?;
    let mut current_record: BamRecord    = BamRecord::default();
    let mut r1_record: Option<BamRecord> = None;
    let mut r2_record: Option<BamRecord> = None;
    let mut prev_qname: Vec<u8>          = Vec::new();
    let mut model_stats: ModelStats      = ModelStats::new(read_length, args.max_ins as f64);

    let process_pair = |r1: &BamRecord, r2: &BamRecord, stats: &mut ModelStats, isizes_modeled: &mut usize| {
        // 1. Accumulate Read Lengths
        stats.r1_lengths.increment(r1.qual.len() as f64);
        stats.r2_lengths.increment(r2.qual.len() as f64);

        // 2. Accumulate cycle-specific base quality distributions (capped at read_length)
        for (cycle, &q) in r1.qual.iter().enumerate().take(read_length) {
            stats.r1_qualities[cycle].increment(q as f64);
        }
        for (cycle, &q) in r2.qual.iter().enumerate().take(read_length) {
            stats.r2_qualities[cycle].increment(q as f64);
        }

        // 3. Filter and accumulate insert sizes
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
                    stats.insert_sizes.increment(insert_size as f64);
                    *isizes_modeled += 1;
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
                    process_pair(r1, r2, &mut model_stats, &mut insert_sizes_modeled);
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

    // Parse the last record pair in the bam file
    if !prev_qname.is_empty() {
        total_pairs += 1;
        if let (Some(r1), Some(r2)) = (&r1_record, &r2_record) {
            process_pair(r1, r2, &mut model_stats, &mut insert_sizes_modeled);
        }
    }

    // Check if processed data
    if total_pairs > 0 && model_stats.r1_lengths.total_count() == 0 {
        return Err(Error::new(
            ErrorKind::InvalidData,
            "No valid pairs found. Is the BAM file name-sorted (`samtools sort -n`)?",
        ));
    }

    // Trim trailing empty bins from all histograms to compress the payload
    model_stats.trim();

    // Write model directly to file via serde_json
    let out_path: &str = &args.model;
    let file: File = File::create(out_path)?;
    serde_json::to_writer_pretty(file, &model_stats).map_err(|e| {
        Error::new(ErrorKind::Other, format!("Failed to save JSON model: {}", e))
    })?;

    // Output a summary to stdout
    let duration: std::time::Duration = start_time.elapsed();
    let summary: serde_json::Value = serde_json::json!({
        "total_pairs_scanned":  total_pairs,
        "insert_sizes_modeled": insert_sizes_modeled,
        "runtime_seconds":      duration.as_secs_f64()
    });

    println!("{}", serde_json::to_string_pretty(&summary).unwrap());

    Ok(())
}