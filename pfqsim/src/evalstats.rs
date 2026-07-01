//! n2bio/pfqsim/src/evalstats.rs
//! 

#![allow(unused)]
use std::collections::HashMap;
use std::io::{self, Error, ErrorKind};

use n2core::bam::{ BamHeader, BamRecord, BamFlags, BamStats };
use crate::cli::AnalyzeArgs;
use crate::config::AnalyzeConfig;

/// Metric storage for an individual alignment profile
#[derive(Default, Debug, Clone)]
pub struct SignalProfile {
    pub mapq: Vec<f64>,
    pub align_score: Vec<f64>,
    pub align_length: Vec<f64>,
    pub as_al: Vec<f64>,
    pub align_proportion: Vec<f64>,
    pub align_accuracy: Vec<f64>,
}

impl SignalProfile {
    /// Inspects a primary record and logs metric points
    pub fn accumulate_record(&mut self, r: &BamRecord) {
        self.mapq.push(r.mapq as f64);

        if let Some(val) = r.get_int_tag(b"AS") {
            self.align_score.push(val as f64);
        }
        if let Some(val) = r.calculate_alignment_length() {
            self.align_length.push(val as f64);
        }
        if let Some(val) = r.calculate_as_al() {
            self.as_al.push(val as f64);
        }
        if let Some(val) = r.calculate_alignment_proportion() {
            self.align_proportion.push(val as f64);
        }
        if let Some(val) = r.calculate_alignment_accuracy() {
            self.align_accuracy.push(val as f64);
        }
    }
}

/// The global tracker holding Positives vs Negatives metrics pools
#[derive(Default, Debug, Clone)]
pub struct GlobalEvaluator {
    /// True Positives: Came from Genome X and safely mapped to Target X
    pub positives: SignalProfile,
    /// False Positives: Came from Genome Y but cross-mapped onto Target X
    pub negatives: SignalProfile,
    
    pub total_expected_simulated: usize,
    pub total_observed_pairs: usize,
}

impl GlobalEvaluator {
    /// Primary business engine: Evaluates a completed pair profile against truth maps
    pub fn evaluate_pair(
        &mut self,
        r1: &BamRecord,
        r2: &BamRecord,
        header: &BamHeader,
        target_to_source: &HashMap<String, String>,
    ) {
        self.total_observed_pairs += 1;

        // 1. Clean the read name: safely slice off the null-terminator (\0) if present
        let mut name_bytes = &r1.read_name[..];
        if let Some(&b'\0') = name_bytes.last() {
            name_bytes = &name_bytes[..name_bytes.len() - 1];
        }
        
        let qname_str = String::from_utf8_lossy(name_bytes);
        
        // Extract true origin from QNAME (Splits "Ecoli_K12:read_42" -> "Ecoli_K12")
        let source_id = match qname_str.split(':').next() {
            Some(id) => id,
            None => return, // Drop malformed query headers safely
        };

        // 2. Validate reference bounds and map the ref_id integer to a string name
        if r1.ref_id < 0 || (r1.ref_id as usize) >= header.references.len() {
            return; // Unmapped (-1) or invalid reference index, drop safely
        }
        let mapped_target_name = &header.references[r1.ref_id as usize].name;

        // 3. Verification Logic
        if let Some(expected_source_for_target) = target_to_source.get(mapped_target_name) {
            if expected_source_for_target == source_id {
                // True Positive Alignment Match
                self.positives.accumulate_record(r1);
                self.positives.accumulate_record(r2);
            } else {
                // False Positive Cross-Contamination Event
                self.negatives.accumulate_record(r1);
                self.negatives.accumulate_record(r2);
            }
        }
    }
}