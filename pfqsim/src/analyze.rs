//! n2bio/pfqsim/src/compose.rs
//! 

use std::io;
use std::fs::File;
use std::path::{ Path, PathBuf };
use std::borrow::Cow;
use std::collections::HashMap;

use n2core::bam::{ BamReader, BamHeader, BamRecord, BamFlags, BamStats };

use crate::config::AnalyzeConfig;
use crate::cli::{ AnalyzeArgs, MappingMode };
use crate::analyze_report::{ AnalyzeReportData, MetricPayload, generate_evaluation_reports };

    
pub(crate) fn run(args: AnalyzeArgs) -> io::Result<()> {
    // Read config
    println!("Loading benchmark baseline configuration: {}", args.config);
    let config: AnalyzeConfig = AnalyzeConfig::from_tsv(&args.config)?;

    println!("Loading reference mapping classifications: {}", args.reference_map);
    let reference_map: HashMap<String, String> = load_reference_map(args.reference_map)?;

    // Total expected targets (excluding negative controls)
    let total_expected_positives: usize = config.rows
        .iter()
        .filter(|row| !row.keyword.eq_ignore_ascii_case("negative"))
        .map(|row| row.calculated_reads)
        .sum();

    // Total expected negative controls
    let total_expected_negatives: usize = config.rows
        .iter()
        .filter(|row| row.keyword.eq_ignore_ascii_case("negative"))
        .map(|row| row.calculated_reads)
        .sum();

    // Initialize the updated Report Data payload
    let mut report_data: AnalyzeReportData = AnalyzeReportData {
        report_name: args.output.clone(),
        total_expected_positives,
        total_expected_negatives,
        total_observed: 0,
        unmapped_true_negatives: 0,
        unmapped_false_negatives: 0,
        true_positives: MetricPayload::new(),
        false_positives_control: MetricPayload::new(),
        false_positives_cross: MetricPayload::new(),
    };

    println!("Processing name-sorted BAM: {}", args.bam);
    let mut bam_reader: BamReader = BamReader::open(&args.bam)?;
    let header: BamHeader = bam_reader.read_header()?;

    let mut current_record: BamRecord = BamRecord::default();
    let mut r1_record: Option<BamRecord> = None;
    let mut r2_record: Option<BamRecord> = None;
    let mut prev_qname: Vec<u8> = Vec::new();

    while bam_reader.read_record(&mut current_record)? {
        if current_record.read_name != prev_qname {
            if !prev_qname.is_empty() {
                if let (Some(r1), Some(r2)) = (&r1_record, &r2_record) {
                    report_data.total_observed += 1;
                    evaluate_and_accumulate_pair(&r1, &r2, &header, &reference_map, &args.mapping_mode, &mut report_data);
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
            report_data.total_observed += 1;
            evaluate_and_accumulate_pair(&r1, &r2, &header, &reference_map, &args.mapping_mode, &mut report_data);
        }
    }

    // Output target creation rules
    let mut out_html_path: PathBuf = Path::new(&args.output).to_path_buf();
    if out_html_path.extension().is_none() {
        out_html_path.set_extension("html");
    }
    
    // Create accompanying JSON file path
    let mut out_json_path: PathBuf = out_html_path.clone();
    out_json_path.set_extension("json");

    println!("Processing metric histograms and writing HTML/JSON reports ...");
    report_data.trim();
    generate_evaluation_reports(&report_data, &out_html_path, &out_json_path)?;
    
    println!("Evaluation completed successfully.\nHTML Dashboard: {}\nRaw JSON Data: {}", 
        out_html_path.display(), 
        out_json_path.display()
    );
    Ok(())
}

/// Loads a 2-column TSV mapping BAM reference targets to their classification group.
fn load_reference_map(path: String) -> io::Result<HashMap<String, String>> {
    let file: File = File::open(PathBuf::from(path))?;
    let mut rdr: csv::Reader<File> = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_reader(file);

    let mut map: HashMap<String, String> = HashMap::new();
    for result in rdr.records() {
        let record: csv::StringRecord = result.map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        if record.len() >= 2 {
            map.insert(record[0].to_string(), record[1].to_string());
        }
    }
    Ok(map)
}

/// Determines if a pair matches the ground-truth assignment and increments scores in the appropriate target histograms
fn evaluate_and_accumulate_pair(
    r1: &BamRecord,
    r2: &BamRecord,
    header: &BamHeader,
    reference_map: &HashMap<String, String>,
    mapping_mode: &MappingMode,
    report_data: &mut AnalyzeReportData,
) {
    // Get read name (query)
    let mut name_bytes: &[u8] = &r1.read_name[..];
    if let Some(&b'\0') = name_bytes.last() {
        name_bytes = &name_bytes[..name_bytes.len() - 1];
    }
    let qname_str: Cow<'_, str> = String::from_utf8_lossy(name_bytes);

    // Parse the read name - "@id:keyword:accession:read_number"
    let mut tokens: std::str::Split<'_, char> = qname_str.split(':');
    let source_id: &str = match tokens.next() {
        Some(val) => val,
        None => return, 
    };
    let keyword: &str = match tokens.next() {
        Some(val) => val,
        None => return, 
    };
    let true_accession: &str = match tokens.next() {
        Some(val) => val,
        None => return, 
    };

    // Determine the expected ground-truth value based on mapping mode
    let expected_truth_value: &str = match mapping_mode {
        MappingMode::ReadId => source_id,
        MappingMode::ReadKeyword => keyword,
        MappingMode::ReadAccession => true_accession,
    };

    let is_negative_control: bool = expected_truth_value.eq_ignore_ascii_case("negative");

    // Check if the read is unmapped (or invalid ref_id)
    if r1.ref_id < 0 || (r1.ref_id as usize) >= header.references.len() {
        if is_negative_control {
            report_data.unmapped_true_negatives += 1;
        } else {
            report_data.unmapped_false_negatives += 1;
        }
        return; // No alignment metrics to capture for unmapped reads
    }

    // Read mapped successfully, evaluate correctness
    let mapped_target_name: &String = &header.references[r1.ref_id as usize].name;

    let target_pool: &mut MetricPayload = if is_negative_control {
        // Negative control that aligned -> FP
        &mut report_data.false_positives_control
    } else if let Some(mapped_classification) = reference_map.get(mapped_target_name) {
        if mapped_classification == expected_truth_value {
            // Correct target mapping -> TP
            &mut report_data.true_positives
        } else {
            // Wrong target mapping -> FP
            &mut report_data.false_positives_cross
        }
    } else {
        // Mapped to a reference not in mapping -> FP
        &mut report_data.false_positives_cross 
    };

    // Increment histogram bins directly via the `.increment()` method
    for record in &[r1, r2] {
        target_pool.mapq.increment(record.mapq as f64);
        
        if let Some(val) = record.get_int_tag(b"AS") {
            target_pool.align_score.increment(val as f64);
        }
        if let Some(val) = record.calculate_alignment_length() {
            target_pool.align_length.increment(val as f64);
        }
        if let Some(val) = record.calculate_as_al() {
            target_pool.as_al.increment(val as f64);
        }
        if let Some(val) = record.calculate_alignment_proportion() {
            target_pool.align_proportion.increment(val as f64);
        }
        if let Some(val) = record.calculate_alignment_accuracy() {
            target_pool.align_accuracy.increment(val as f64);
        }
    }
}