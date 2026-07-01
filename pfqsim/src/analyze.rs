//! n2bio/pfqsim/src/compose.rs
//! 

use std::io;
use std::path::{ Path, PathBuf };
use std::borrow::Cow;

use n2core::bam::{ BamReader, BamHeader, BamRecord, BamFlags, BamStats };

use crate::config::AnalyzeConfig;
use crate::cli::AnalyzeArgs;
use crate::report::{ AnalyzeReportData, MetricPayload, generate_evaluation_report };

    
pub(crate) fn run(args: AnalyzeArgs) -> io::Result<()> {
    // Read config
    println!("Loading benchmark baseline configuration: {}", args.manifest);
    let config: AnalyzeConfig = AnalyzeConfig::from_tsv(&args.manifest)?;

    // Total simulated inputs to calibrate baseline False Negatives
    let total_expected: usize = config.rows.iter().map(|row| row.reads).sum();

    let mut report_data: AnalyzeReportData = AnalyzeReportData {
        report_name: args.output.clone(),
        total_expected_pairs: total_expected,
        total_observed_pairs: 0,
        positives: MetricPayload::default(),
        negatives: MetricPayload::default(),
    };

    println!("Processing name-sorted BAM tracking pipeline: {}", args.bam.display());
    let mut bam_reader: BamReader = BamReader::open(args.bam.to_str().unwrap())?;
    let header: BamHeader = bam_reader.read_header()?;

    let mut current_record: BamRecord = BamRecord::default();
    let mut r1_record: Option<BamRecord> = None;
    let mut r2_record: Option<BamRecord> = None;
    let mut prev_qname: Vec<u8> = Vec::new();

    while bam_reader.read_record(&mut current_record)? {
        if current_record.read_name != prev_qname {
            if !prev_qname.is_empty() {
                if let (Some(r1), Some(r2)) = (&r1_record, &r2_record) {
                    report_data.total_observed_pairs += 1;
                    evaluate_and_accumulate_pair(&r1, &r2, &header, &mut report_data);
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

    // Process final records
    if !prev_qname.is_empty() {
        if let (Some(r1), Some(r2)) = (&r1_record, &r2_record) {
            report_data.total_observed_pairs += 1;
            evaluate_and_accumulate_pair(&r1, &r2, &header, &mut report_data);
        }
    }

    // Output target creation rules
    let mut out_html_path: PathBuf = Path::new(&args.output).to_path_buf();
    if out_html_path.extension().is_none() {
        out_html_path.set_extension("html");
    }

    println!("Processing metric data vectors and writing HTML Report dashboard...");
    generate_evaluation_report(&mut report_data, &out_html_path)?;
    
    println!("Evaluation completed successfully. Report dashboard ready at: {}", out_html_path.display());
    Ok(())
}

/// Determines if a pair matches the ground-truth assignment and moves scores into the appropriate target vectors
fn evaluate_and_accumulate_pair(
    r1: &BamRecord,
    r2: &BamRecord,
    header: &BamHeader,
    report_data: &mut AnalyzeReportData,
) {
    // Get read name (query)
    let mut name_bytes: &[u8] = &r1.read_name[..];
    if let Some(&b'\0') = name_bytes.last() {
        name_bytes = &name_bytes[..name_bytes.len() - 1];
    }
    let qname_str: Cow<'_, str> = String::from_utf8_lossy(name_bytes);

    // Parse the read name - "@id:accession:read_number"
    let mut tokens: std::str::Split<'_, char> = qname_str.split(':');
    let _source_id: &str = match tokens.next() {
        Some(id) => id,
        None => return, // Malformed name string, skip
    };
    let true_accession: &str = match tokens.next() {
        Some(acc) => acc,
        None => return, // Malformed name string, skip
    };

    // Get the alignment target name
    if r1.ref_id < 0 || (r1.ref_id as usize) >= header.references.len() {
        return; 
    }
    let mapped_target_accession: &String = &header.references[r1.ref_id as usize].name;

    // Evaluate alignment
    let target_pool: &mut MetricPayload = if true_accession == mapped_target_accession {
        &mut report_data.positives
    } else {
        &mut report_data.negatives
    };

    // Append alignment data profiles into vectors
    for record in &[r1, r2] {
        target_pool.mapq.push(record.mapq as f64);
        
        if let Some(val) = record.get_int_tag(b"AS") {
            target_pool.align_score.push(val as f64);
        }
        if let Some(val) = record.calculate_alignment_length() {
            target_pool.align_length.push(val as f64);
        }
        if let Some(val) = record.calculate_as_al() {
            target_pool.as_al.push(val as f64);
        }
        if let Some(val) = record.calculate_alignment_proportion() {
            target_pool.align_proportion.push(val as f64);
        }
        if let Some(val) = record.calculate_alignment_accuracy() {
            target_pool.align_accuracy.push(val as f64);
        }
    }
}