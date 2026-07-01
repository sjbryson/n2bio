//! n2bio/pfqsim/src/compose.rs
//! 

#![allow(unused)]

use std::io::{self, BufRead, Write, Error, ErrorKind};
use std::collections::HashMap;
use std::fs::File;

use crate::config::AnalyzeConfig;
use crate::cli::AnalyzeArgs;
use crate::evalstats::GlobalEvaluator;
use n2core::bam::{ BamReader, BamHeader, BamRecord, BamFlags };

    
pub(crate) fn run(args: AnalyzeArgs) -> io::Result<()> {
    // 1. Read config
    println!("Loading analysis baseline configurations: {}", args.manifest);
    let config: AnalyzeConfig = AnalyzeConfig::from_tsv(&args.manifest)?;

    // 2. Set up read-to-target map and get total read counts
    let mut target_to_source: HashMap<String, String> = HashMap::new();
    let mut total_expected: usize = 0;

    for row in &config.rows {
        target_to_source.insert(row.target.clone(), row.id.clone());
        total_expected += row.reads;
    }

    let mut evaluator: GlobalEvaluator = GlobalEvaluator {
        total_expected_simulated: total_expected,
        ..Default::default()
    };

    // 3. Initiate BAM Parser Loop
    println!("Evaluating BAM alignment pairs: {}", args.bam.display());
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
                    evaluator.evaluate_pair(&r1, &r2, &header, &target_to_source);
                }
            }
            r1_record = None;
            r2_record = None;
            prev_qname = current_record.read_name.clone();
        }

        // Just use primary alignments 
        if current_record.is_primary() { //&& current_record.is_proper() { 
            if current_record.is_read1() {
                r1_record = Some(current_record.clone());
            } else if current_record.is_read2() {
                r2_record = Some(current_record.clone());
            }
        }
    }

    // Capture the trailing record pair remaining inside buffer streams
    if !prev_qname.is_empty() {
        if let (Some(r1), Some(r2)) = (&r1_record, &r2_record) {
            evaluator.evaluate_pair(&r1, &r2, &header, &target_to_source);
        }
    }

    // Print summary snapshot of metrics gathered so far
    println!("------------------------------------------------------------");
    println!("Analysis Complete!");
    println!("Simulated Ground Truth Pairs : {}", evaluator.total_expected_simulated);
    println!("Observed BAM Primary Pairs   : {}", evaluator.total_observed_pairs);
    println!("True Positive Data points    : {}", evaluator.positives.mapq.len());
    println!("False Positive Data points   : {}", evaluator.negatives.mapq.len());
    println!("------------------------------------------------------------");

    Ok(())
}