//! n2bio/peat/src/binreads.rs
//! 

use std::io;
use n2bio::bam::{ BamReader, BamHeader, BamRecord, BamFlags };
use n2bio::fastq::{ PairedFastqRecord, PairedRead};

use crate::cli::{ BinReadsArgs, ThresholdMetrics };
use crate::alignmentfilters::{ highpass_bamfilter, threshold_args };
use crate::binreads_pooler::BinnedFastqPool;
use crate::binreads_resolver::BinResolver;
use crate::binreads_stats::BinReadReport;


pub(crate) fn run(args: BinReadsArgs) -> io::Result<()> {
    // Initialize BamReader and BamHeader
    let mut bam_reader: BamReader = BamReader::open(args.bam.to_str().unwrap())?;
    let header: BamHeader = bam_reader.read_header()?;

    // Initialize target resolver
    let bin_resolver: BinResolver = BinResolver::new(args.reference_map.as_deref())?;

    // Setup threshold evaluator state
    let threshold_active: bool = threshold_args(&args.thresholds);

    // Setup Writer Pool
    let max_open_bins: usize = 64;
    let gz_threads: usize = if args.threads > 2 { 2 } else { 1 };

    let mut pool: BinnedFastqPool = BinnedFastqPool::new(
        max_open_bins,
        1_000, // Buffer 1,000 paired reads per bin before disk flush
        args.output_dir.clone(),
        gz_threads,
        args.report.clone(),
    );

    let mut bin_stats: BinReadReport = BinReadReport::default();

    // Streaming Loop State
    let mut current_record: BamRecord = BamRecord::default();
    let mut r1_record: Option<BamRecord> = None;
    let mut r2_record: Option<BamRecord> = None;
    let mut prev_qname: Vec<u8> = Vec::new();

    // Main Streaming Reader Loop (Name-Sorted Expectation)
    while bam_reader.read_record(&mut current_record)? {
        // When qname changes, evaluate the previous read pair
        if current_record.read_name != prev_qname {
            if let (Some(r1), Some(r2)) = (r1_record.take(), r2_record.take()) {
                evaluate_and_bin_pair(
                    r1,
                    r2,
                    &header,
                    &bin_resolver,
                    &args.thresholds,
                    threshold_active,
                    &mut pool,
                    &mut bin_stats,
                )?;
            }
            prev_qname = current_record.read_name.clone();
        }

        // Only buffer primary alignments for pairing
        if current_record.is_primary() {
            if current_record.is_read1() {
                r1_record = Some(current_record.clone());
            } else if current_record.is_read2() {
                r2_record = Some(current_record.clone());
            }
        }
    }

    // Flush trailing final pair after EOF
    if let (Some(r1), Some(r2)) = (r1_record.take(), r2_record.take()) {
        evaluate_and_bin_pair(
            r1,
            r2,
            &header,
            &bin_resolver,
            &args.thresholds,
            threshold_active,
            &mut pool,
            &mut bin_stats,
        )?;
    }

    // Flush all writer buffers and close gzip handles
    pool.finish_all()?;

    // Write JSON report
    let report_path: std::path::PathBuf = args.output_dir.join(format!("{}.bin_report.json", args.report));
    bin_stats.write_json(&report_path)?;

    Ok(())
}


/// Helper function to filter, orient, resolve bins, and write a pair to disk.
fn evaluate_and_bin_pair(
    r1: BamRecord,
    r2: BamRecord,
    header: &BamHeader,
    resolver: &BinResolver,
    thresholds: &ThresholdMetrics,
    threshold_active: bool,
    pool: &mut BinnedFastqPool,
    stats: &mut BinReadReport,
) -> io::Result<()> {
    // Highpass filter
    let r1_pass: bool = highpass_bamfilter(&r1, thresholds, threshold_active);
    let r2_pass: bool = highpass_bamfilter(&r2, thresholds, threshold_active);

    if !r1_pass && !r2_pass {
        return Ok(());
    }

    // Resolve Target Names -> Bin IDs
    let bin_r1: Option<&str> = if r1_pass {
        r1.target_name(header).and_then(|t| resolver.get_bin(t))
    } else {
        None
    };

    let bin_r2: Option<&str> = if r2_pass {
        r2.target_name(header).and_then(|t| resolver.get_bin(t))
    } else {
        None
    };

    // Convert BamRecord -> PairedRead -> PairedFastqRecord
    let p_r1: PairedRead = PairedRead::from_bamrec(&r1);
    let p_r2: PairedRead = PairedRead::from_bamrec(&r2);

    let (rec_r1, rec_r2) = match (p_r1, p_r2) {
        (PairedRead::R1(a), PairedRead::R2(b)) => (a, b),
        (PairedRead::R2(b), PairedRead::R1(a)) => (a, b),
        _ => return Ok(()), // Safeguard against malformed flags
    };

    let pair: PairedFastqRecord = PairedFastqRecord { r1: rec_r1, r2: rec_r2 };

    // Route pair based on bin assignments
    match (bin_r1, bin_r2) {
        (Some(b1), Some(b2)) if b1 == b2 => {
            stats.inc_concordant(b1);
            pool.push(b1, pair)?;
        }
        (Some(b1), Some(b2)) => {
            stats.inc_discordant(b1);
            stats.inc_discordant(b2);
            pool.push(b1, pair)?; // Primary bin write
        }
        (Some(b1), None) => {
            stats.inc_r1_orphan(b1);
            pool.push(b1, pair)?;
        }
        (None, Some(b2)) => {
            stats.inc_r2_orphan(b2);
            pool.push(b2, pair)?;
        }
        (None, None) => {}
    }

    Ok(())
}