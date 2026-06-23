//! n2core/bamrep/src/main.rs
//! 

use clap::Parser;
use std::collections::HashMap;
use std::fs::File;
use std::path::PathBuf;
use rayon::prelude::*;

mod stats;
mod report;
use n2core::bam::{ BamReader, BamHeader, BamRecord, BamFlags };
use crate::stats::{ GlobalStats, StatSummary, StatsAccumulator, ReportData };

// ============================================================================
// Args
// ============================================================================

#[derive(Parser, Debug)]
#[command(author, version, about = "BAM Alignment Stats")]
struct Args {
    /// Input name-sorted BAM file
    #[arg(short = 'b', long)]
    bam: PathBuf,

    /// Report file prefix - creates {report}.json and optional {report}.html
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
    // Variables for primary alignment pairing
    let mut current_record: BamRecord = BamRecord::default();
    let mut r1_record: Option<BamRecord> = None;
    let mut r2_record: Option<BamRecord> = None;
    let mut prev_qname: Vec<u8> = Vec::new();
    let mut total_pairs: i32 = 0;
    let mut global_stats: GlobalStats = GlobalStats::default();
    let mut alignment_stats: StatsAccumulator = StatsAccumulator::default();
   
    // Read loop
    while bam_reader.read_record(&mut current_record)? {
        global_stats.update_counts(&current_record);
        alignment_stats.update_read_stats(&current_record, args.max_len);

        if current_record.read_name != prev_qname {
            if !prev_qname.is_empty() {
                total_pairs += 1;
                // Only calculate insert size if valid primary/proper pair for this name
                if let (Some(r1), Some(r2)) = (&r1_record, &r2_record) {
                    alignment_stats.update_insert_size(r1, r2, args.min_mapq, args.max_ins as i32);
                }
            }
            r1_record = None;
            r2_record = None;
            prev_qname = current_record.read_name.clone();
        }

        if current_record.is_primary() && current_record.is_proper() {
            if current_record.is_read1() {
                r1_record = Some(current_record.clone());
            } else if current_record.is_read2() {
                r2_record = Some(current_record.clone());
            }
        }

    }

    // Process final record pair
    if !prev_qname.is_empty() {
        total_pairs += 1;
        if let (Some(r1), Some(r2)) = (&r1_record, &r2_record) {
            alignment_stats.update_insert_size(r1, r2, args.min_mapq, args.max_ins as i32);
        }
    }

    if total_pairs == 0 {
        eprintln!("No valid pairs found. Is the BAM file name-sorted (`samtools sort -n`)?");
        std::process::exit(1);
    }

    //Update GlobalStats with the total_pairs count
    global_stats.set_total_reads(total_pairs as usize);

    println!("BAM reading complete. Processed {} pairs. Generating summaries and plots...", total_pairs);

    // Group all 11 vectors into a list so Rayon can process them concurrently
    let mut stats_to_process: Vec<(&str, &mut Vec<f64>)> = vec![
        ("pe_insert_size", &mut alignment_stats.pe_insert_size),
        ("r1_mapq", &mut alignment_stats.r1_mapq),
        ("r2_mapq", &mut alignment_stats.r2_mapq),
        ("r1_align_score", &mut alignment_stats.r1_align_score),
        ("r2_align_score", &mut alignment_stats.r2_align_score),
        ("r1_align_length", &mut alignment_stats.r1_align_length),
        ("r2_align_length", &mut alignment_stats.r2_align_length),
        ("r1_as_al", &mut alignment_stats.r1_as_al),
        ("r2_as_al", &mut alignment_stats.r2_as_al),
        ("r1_align_proportion", &mut alignment_stats.r1_align_proportion),
        ("r2_align_proportion", &mut alignment_stats.r2_align_proportion),
        ("r1_align_accuracy", &mut alignment_stats.r1_align_accuracy),
        ("r2_align_accuracy", &mut alignment_stats.r2_align_accuracy),
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
    let mut alignment_results: HashMap<String, StatSummary> = HashMap::new();
    for (name, summary) in results_vec {
        alignment_results.insert(name, summary);
    }

    let report_data: ReportData = ReportData {
        global_stats,
        alignment_stats: alignment_results,
    };

    // Write JSON Report
    let mut json_path: PathBuf = args.report.clone();
    json_path.set_extension("json"); // Enforces the .json extension

    let report_file: File = File::create(&json_path)?;
    serde_json::to_writer_pretty(report_file, &report_data)?;

    // Write optional HTML report
    if args.html {
        let mut html_path: PathBuf = args.report.clone();
        html_path.set_extension("html");

        if let Err(e) = report::generate_html_report(&report_data, &html_path, &args.report) {
            eprintln!("Warning: Failed to generate HTML report: {}", e);
        }
    }

    println!("Analysis complete.");
    println!("Processed {} pairs.", total_pairs);
    println!("Results saved to {:?}", &json_path);
    
    Ok(())
}

