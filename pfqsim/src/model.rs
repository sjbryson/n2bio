//! n2bio/pfqsim/src/model.rs
//! 

use std::io::{self, Error, ErrorKind};
use std::time::Instant;
use std::collections::BTreeMap;

use n2core::bam::{ BamReader, BamRecord, CigarOp };

use crate::cli::ModelArgs;
use crate::simstats::{
    LibraryModel, InsertModel, QualityModel, NormalDistParams, update_qscore_model
};

// ============================================================================
// Helper Functions
// ============================================================================

/// Calculates how many reference bases the alignment spans using the CIGAR string
fn calculate_ref_span(cigar: &[(CigarOp, u32)]) -> i32 {
    cigar.iter()
        .filter(|(op, _)| matches!(
            op,
            CigarOp::Match | CigarOp::Deletion | CigarOp::Skip | 
            CigarOp::SequenceMatch | CigarOp::SequenceMismatch
        ))
        .map(|(_, len)| *len as i32)
        .sum()
}

// ============================================================================
// Main Runner
// ============================================================================

pub fn run(args: ModelArgs) -> io::Result<()> {
    let start_time: Instant             = Instant::now();
    let read_length: usize              = args.length;
    let mut total_pairs: usize          = 0usize;
    let mut insert_sizes_modeled: usize = 0usize;
    
    // Use temporary BTreeMaps to accumulate the frequencies.
    let mut raw_insert_sizes: BTreeMap<i32, usize> = BTreeMap::new();
    let mut raw_r1_qual: Vec<BTreeMap<u8, usize>>  = Vec::new();
    let mut raw_r2_qual: Vec<BTreeMap<u8, usize>>  = Vec::new();
    
    let mut bam_reader: BamReader        = BamReader::open(args.bam.to_str().unwrap())?;
    let _header: n2core::bam::BamHeader  = bam_reader.read_header()?;
    let mut current_record: BamRecord    = BamRecord::default();
    let mut r1_record: Option<BamRecord> = None;
    let mut r2_record: Option<BamRecord> = None;
    let mut prev_qname: Vec<u8>          = Vec::new();

    // A closure that mutates our raw accumulators
    let mut process_pair = |r1: &BamRecord, r2: &BamRecord, sizes_count: &mut usize| {
        update_qscore_model(&mut raw_r1_qual, &r1.qual);
        update_qscore_model(&mut raw_r2_qual, &r2.qual);

        if r1.mapq as usize >= args.mapq && r2.mapq as usize >= args.mapq {
            if r1.ref_id == r2.ref_id && r1.ref_id != -1 {
                let (fwd, rev) = if r1.pos <= r2.pos {
                    (r1, r2)
                } else {
                    (r2, r1)
                };

                let rev_cigar: Vec<(CigarOp, u32)> = rev.parsed_cigar();
                let ref_span: i32 = calculate_ref_span(&rev_cigar);
                let insert_size: i32 = (rev.pos + ref_span) - fwd.pos;
                
                if insert_size > 0 {
                    *raw_insert_sizes.entry(insert_size).or_insert(0) += 1;
                    *sizes_count += 1;
                }
            }
        }
    };

    while bam_reader.read_record(&mut current_record)? {
        if current_record.read_name != prev_qname {
            if !prev_qname.is_empty() {
                total_pairs += 1;
                if let (Some(r1), Some(r2)) = (&r1_record, &r2_record) {
                    process_pair(r1, r2, &mut insert_sizes_modeled);
                }
            }
            r1_record  = None;
            r2_record  = None;
            prev_qname = current_record.read_name.clone();
        }

        if current_record.flag & 0x40 != 0 {
            r1_record = Some(current_record.clone());
        } else if current_record.flag & 0x80 != 0 {
            r2_record = Some(current_record.clone());
        }
    }

    // Parse last record in the file
    if !prev_qname.is_empty() {
        total_pairs += 1;
        if let (Some(r1), Some(r2)) = (&r1_record, &r2_record) {
            process_pair(r1, r2, &mut insert_sizes_modeled);
        }
    }

    if total_pairs > 0 && raw_r1_qual.is_empty() {
        return Err(Error::new(
            ErrorKind::InvalidData,
            "No valid pairs found. Is the BAM file name-sorted (`samtools sort -n`)?",
        ));
    }

    // Convert raw Q-score counts into final NormalDistParams per position
    let process_qscores = |raw_q: &[BTreeMap<u8, usize>]| -> Vec<NormalDistParams> {
        let mut params: Vec<NormalDistParams> = Vec::with_capacity(read_length);
        for i in 0..read_length {
            if i < raw_q.len() {
                params.push(NormalDistParams::from_frequency_map(&raw_q[i]));
            } else {
                let last_known: NormalDistParams = params.last().cloned().unwrap_or(NormalDistParams { mean: 30.0, std_dev: 2.0 });
                params.push(last_known);
            }
        }
        params
    };

    // Assemble the final nested model structure
    let final_model: LibraryModel = LibraryModel {
        insert_size: InsertModel {
            insert_dist: NormalDistParams::from_frequency_map(&raw_insert_sizes),
        },
        quality: QualityModel {
            r1_quals: process_qscores(&raw_r1_qual),
            r2_quals: process_qscores(&raw_r2_qual),
        },
    };

    // Write models .to_file()
    let out_path = args.output.to_str().ok_or_else(|| {
        Error::new(ErrorKind::InvalidInput, "Output path is not valid UTF-8")
    })?;
    
    final_model.to_file(out_path).map_err(|e| {
        Error::new(ErrorKind::Other, format!("Failed to save JSON model: {}", e))
    })?;

    let duration: std::time::Duration = start_time.elapsed();
    let summary: serde_json::Value = serde_json::json!({
        "total_pairs_scanned": total_pairs,
        "insert_sizes_modeled": insert_sizes_modeled,
        "runtime_seconds": duration.as_secs_f64()
    });

    println!("{}", serde_json::to_string_pretty(&summary).unwrap());

    Ok(())
}