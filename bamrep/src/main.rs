//! n2core/bamrep/src/main.rs
//! 

use clap::Parser;
use std::collections::HashMap;
use std::fs::File;
use std::path::PathBuf;
use rayon::prelude::*;

mod stats;
mod report;
use n2core::bam::{ BamReader, BamHeader, BamRecord, BamStats };
use crate::stats::{ StatSummary, StatsAccumulator };

// ============================================================================
// Args
// ============================================================================

#[derive(Parser, Debug)]
#[command(author, version, about = "BAM Alignment Stats")]
struct Args {
    /// Input name-sorted BAM file
    #[arg(short = 'b', long)]
    bam: PathBuf,

    /// Output JSON report file
    #[arg(short = 'r', long)]
    report: PathBuf,

    /// Generate html plots
    #[arg(long)]
    html: bool,

    /// Minimum MAPQ score for insert size calculation
    #[arg(short = 'q', long, default_value_t = 40)]
    min_mapq: usize,

    /// Max insert size to use for summary stats calculation
    #[arg(short = 'i', long, default_value_t = 1000)]
    max_ins: usize,

    /// Max read length to use
    #[arg(short = 'l', long, default_value_t = 150)]
    max_len: usize,
}

// ============================================================================
// Main
// ============================================================================

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Args = Args::parse();
    // Setup BamReader
    let mut bam_reader: BamReader = BamReader::open(args.bam.to_str().unwrap())?;
    let _header: BamHeader = bam_reader.read_header()?;
    // Variables for BAM parsing
    let mut current_record: BamRecord = BamRecord::default();
    let mut r1_record: Option<BamRecord> = None;
    let mut r2_record: Option<BamRecord> = None;
    let mut prev_qname: Vec<u8> = Vec::new();
    // Global stats variables
    let mut total_pairs: i32 = 0;
    let mut stats: StatsAccumulator = StatsAccumulator::default();
    // Set some max parameter values based on args.max_len
    let max_mapq: f64 = 60.0;
    let max_as: f64 = 2.0 * args.max_len as f64;
    let max_al: f64 = args.max_len as f64;
    let max_as_al: f64 = 2.0;
    let max_align_prop: f64 = 1.0;
    let max_align_pi: f64 = 100.0;
    // Closure to extract alignment stats from individual reads
    let mut extract_read_stats = |r: &BamRecord, is_r1: bool| {
        // Route the data to the correct vectors
        let (
            mapq_vec, 
            align_score_vec,
            align_len_vec, 
            as_al_vec, 
            align_prop_vec, 
            align_acc_vec,
        ) = if is_r1 {
            (
                &mut stats.r1_mapq,
                &mut stats.r1_align_score, 
                &mut stats.r1_align_length, 
                &mut stats.r1_as_al, 
                &mut stats.r1_align_proportion, 
                &mut stats.r1_align_accuracy
            )
        } else {
            (
                &mut stats.r2_mapq,
                &mut stats.r2_align_score, 
                &mut stats.r2_align_length, 
                &mut stats.r2_as_al, 
                &mut stats.r2_align_proportion, 
                &mut stats.r2_align_accuracy
            )
        };  
        // If read is unmapped push 0 to each vec
        //if r.mapq == 0 {
        //    mapq_vec.push(0 as f64);
        //    align_score_vec.push(0 as f64);
        //    align_len_vec.push(0 as f64);
        //    as_al_vec.push(0 as f64);
        //    align_prop_vec.push(0 as f64);
        //    align_acc_vec.push(0 as f64);
        //}

        // Only push alignment metrics if the read is mapped
        if r.mapq > 0 {
            mapq_vec.push((r.mapq as f64).min(max_mapq));
            
            if let Some(val) = r.get_int_tag(b"AS") { 
                align_score_vec.push((val as f64).min(max_as)); 
            }
            if let Some(val) = r.calculate_alignment_length() { 
                align_len_vec.push((val as f64).min(max_al)); 
            }
            if let Some(val) = r.calculate_as_al() { 
                as_al_vec.push((val as f64).min(max_as_al)); 
            }
            if let Some(val) = r.calculate_alignment_proportion() { 
                align_prop_vec.push((val as f64).min(max_align_prop)); 
            }
            if let Some(val) = r.calculate_alignment_accuracy() { 
                align_acc_vec.push((val as f64).min(max_align_pi)); 
            }
        }
    };

    // Closure to process the full pair
    let mut process_pair = |r1: &BamRecord, r2: &BamRecord| {
        extract_read_stats(r1, true);
        extract_read_stats(r2, false);

        // Mapq filtering and insert size logic - only use good alignments
        if r1.mapq as usize >= args.min_mapq && r2.mapq as usize >= args.min_mapq {
            if r1.ref_id == r2.ref_id && r1.ref_id != -1 {
                let (fwd, rev) = if r1.pos <= r2.pos { (r1, r2) } else { (r2, r1) };
                let ref_span: i32 = rev.calculate_ref_span().unwrap_or(0) as i32;
                let insert_size: i32 = (rev.pos + ref_span) - fwd.pos;
                
                if insert_size > 0 && insert_size <= args.max_ins as i32 {
                    stats.pe_insert_size.push(insert_size as f64);
                }
            }
        }
    };

    // Read loop
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

        // FLAG 0x40 = first in pair (R1), FLAG 0x80 = second in pair (R2)
        if current_record.flag & 0x40 != 0 {
            r1_record = Some(current_record.clone());
        } else if current_record.flag & 0x80 != 0 {
            r2_record = Some(current_record.clone());
        }
    }

    // Process final record pair
    if !prev_qname.is_empty() {
        total_pairs += 1;
        if let (Some(r1), Some(r2)) = (&r1_record, &r2_record) {
            process_pair(r1, r2);
        }
    }

    if total_pairs == 0 {
        eprintln!("No valid pairs found. Is the BAM file name-sorted (`samtools sort -n`)?");
        std::process::exit(1);
    }

   println!("BAM reading complete. Processed {} pairs. Generating summaries and plots...", total_pairs);

    // Group all 11 vectors into a list so Rayon can process them concurrently
    let mut stats_to_process: Vec<(&str, &mut Vec<f64>)> = vec![
        ("pe_insert_size", &mut stats.pe_insert_size),
        ("r1_mapq", &mut stats.r1_mapq),
        ("r2_mapq", &mut stats.r2_mapq),
        ("r1_align_score", &mut stats.r1_align_score),
        ("r2_align_score", &mut stats.r2_align_score),
        ("r1_align_length", &mut stats.r1_align_length),
        ("r2_align_length", &mut stats.r2_align_length),
        ("r1_as_al", &mut stats.r1_as_al),
        ("r2_as_al", &mut stats.r2_as_al),
        ("r1_align_proportion", &mut stats.r1_align_proportion),
        ("r2_align_proportion", &mut stats.r2_align_proportion),
        ("r1_align_accuracy", &mut stats.r1_align_accuracy),
        ("r2_align_accuracy", &mut stats.r2_align_accuracy),
    ];

    let max_ins_val: f64 = args.max_ins as f64;
    let max_len_val: f64 = args.max_len as f64;
    // Process all 11 stats concurrently using Rayon
    let results_vec: Vec<(String, StatSummary)> = stats_to_process
        .par_iter_mut()
        .map(|(name, data)| {
            // Pass the name to get the correct config
            let summary: StatSummary = StatSummary::calculate(name, data, max_ins_val, max_len_val);
            (name.to_string(), summary)
        })
        .collect();

    // Reconstruct the HashMap for JSON serialization
    let mut results: HashMap<String, StatSummary> = HashMap::new();
    for (name, summary) in results_vec {
        results.insert(name, summary);
    }

    // Write JSON Report
    let mut json_path: PathBuf = args.report.clone();
    json_path.set_extension("json"); // Enforces the .json extension

    let report_file: File = File::create(&json_path)?;
    serde_json::to_writer_pretty(report_file, &results)?;

    // Write optional HTML report
    if args.html {
        let mut html_path: PathBuf = args.report.clone();
        html_path.set_extension("html");

        if let Err(e) = report::generate_html_report(&results, &html_path) {
            eprintln!("Warning: Failed to generate HTML report: {}", e);
        }
    }


    
    println!("Analysis complete.");
    println!("Processed {} pairs.", total_pairs);
    println!("Results saved to {:?}", &json_path);
    
    Ok(())
}

