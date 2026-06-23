//! n2core/bamrep/src/stats.rs
//! 
use serde::Serialize;
use std::collections::HashMap;

use n2core::bam::{ BamRecord, BamFlags, BamStats };
use crate::report::ReportConfig;

// ============================================================================
// Global Stats
// ============================================================================

// Holds the counts for a single read direction (R1 or R2)
#[derive(Serialize, Default, Clone, Debug)]
pub struct ReadStats {
    pub total_reads: usize,        // Total count of each R1 and R2 read
    pub total_unaligned: usize,    // Total R1 or R2 unaligned reads
    pub total_alignments: usize,   // Total R1 or R2 alignments in BAM file
    pub primary_mapped: usize,     // Total R1 or R2 primary alignments
    pub primary_mapq: usize,       // Total R1 or R2 primary alignments with MAPQ > 0
    pub concordant_mapped: usize,  // Total R1 or R2 concordant primary alignments
    pub discordant_mapped: usize,  // Total R1 or R2 discordant primary alignments
    pub singletons: usize,         // Total R1 or R2 singleton primary alignments
    pub secondary_mapped: usize,   // Total R1 or R2 secondary or accessory alignments
    pub mapq_0: usize,             // Total R1 or R2 alignments with MAPQ = 0
}

// Holds the combined R1/R2 global stats
#[derive(Serialize, Default, Clone, Debug)]
pub struct GlobalStats {
    pub r1: ReadStats,
    pub r2: ReadStats,
}

impl GlobalStats {
    pub fn update_counts(&mut self, r: &BamRecord) {
        // Route to the correct struct r1 or r2
        let mate: &mut ReadStats = if r.is_read1() {
            &mut self.r1
        } else if r.is_read2() {
            &mut self.r2
        } else {
            return; // Skip if neither R1 nor R2 (e.g., unpaired)
        };
        if !r.is_mapped() { mate.total_unaligned += 1; }
        if r.is_mapped() {
            mate.total_alignments += 1;
            if r.mapq == 0 {
                mate.mapq_0 += 1;
            }
            if r.is_primary() {
                mate.primary_mapped += 1;
                if r.mapq > 0 {
                    mate.primary_mapq += 1;
                }
                
                if r.is_mate_unmapped() {
                    mate.singletons += 1;
                } else if r.is_proper() {
                    mate.concordant_mapped += 1;
                } else {
                    mate.discordant_mapped += 1;
                }
            } else {
                mate.secondary_mapped += 1;
            }
        }
    }

    pub fn set_total_reads(&mut self, n: usize) {
        self.r1.total_reads = n.clone();
        self.r2.total_reads = n.clone();
    }
}

// ============================================================================
// Summary Stats
// ============================================================================

// The new root JSON structure
#[derive(Serialize)]
pub struct ReportData {
    pub global_stats: GlobalStats,
    pub alignment_stats: HashMap<String, StatSummary>,
}

#[derive(Serialize, Default)]
pub struct HistogramData {
    bin_min: f64,
    bin_size: f64,
    counts: Vec<u32>,
}

#[derive(Serialize, Default)]
pub struct StatSummary {
    pub count: f64,
    pub mean: f64,
    pub median: f64,
    pub stdev: f64,
    pub min: f64,
    pub max: f64,
    pub histogram: HistogramData,
}

impl StatSummary {
    pub fn calculate(name: &str, data: &mut [f64], max_ins: f64, max_len: f64) -> Self {
        if data.is_empty() {
            return Self::default();
        }
        // Sort data for median, min, and max calculations
        data.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        let count: f64 = data.len() as f64;
        let sum:   f64 = data.iter().sum();
        let mean:  f64 = sum / count;
        let min:   f64 = data[0];
        let max:   f64 = data[data.len() - 1];
        
        let median: f64 = if data.len() % 2 == 0 {
            let mid: usize = data.len() / 2;
            (data[mid - 1] + data[mid]) / 2.0
        } else {
            data[data.len() / 2]
        };

        // population vaiance
        let variance: f64 = data.iter().map(|&value| {
            let diff: f64 = mean - value;
            diff * diff
        }).sum::<f64>() / count;

        let stdev: f64 = variance.sqrt();

        let config: ReportConfig = ReportConfig::from_stat(name, min, max, max_ins, max_len);
        let num_bins: usize = ((config.max - config.min) / config.bin_size).ceil() as usize;
        let mut bins: Vec<u32> = vec![0u32; num_bins.max(1)];

        for &val in data.iter() {
            if val < config.min || val >= config.max { continue; }
            let mut idx: usize = ((val - config.min) / config.bin_size) as usize;
            if idx >= bins.len() { idx = bins.len() - 1; }
            bins[idx] += 1;
        }

        StatSummary {
            count, mean, median, stdev, min, max,
            histogram: HistogramData {
                bin_min: config.min,
                bin_size: config.bin_size,
                counts: bins,
            }
        }
    }
}

/// Holds raw extracted values to compute summaries and plot histograms later
#[derive(Default)]
pub struct StatsAccumulator {
    pub pe_insert_size: Vec<f64>,

    pub r1_mapq: Vec<f64>,
    pub r1_align_score: Vec<f64>,
    pub r1_align_length: Vec<f64>,
    pub r1_as_al: Vec<f64>,
    pub r1_align_proportion: Vec<f64>,
    pub r1_align_accuracy: Vec<f64>,

    pub r2_mapq: Vec<f64>,
    pub r2_align_score: Vec<f64>,
    pub r2_align_length: Vec<f64>,
    pub r2_as_al: Vec<f64>,
    pub r2_align_proportion: Vec<f64>,
    pub r2_align_accuracy: Vec<f64>,
}

impl StatsAccumulator {
    pub fn update_read_stats(&mut self, r: &BamRecord, max_len: usize) {
        // Only process mapped, primary reads with MAPQ > 0
        if r.mapq == 0 || !r.is_primary() || !r.is_mapped() {
            return;
        }

        // Setup max thresholds
        let max_mapq: f64 = 60.0;
        let max_as: f64 = 2.0 * max_len as f64;
        let max_al: f64 = max_len as f64;
        let max_as_al: f64 = 2.0;
        let max_align_prop: f64 = 1.0;
        let max_align_pi: f64 = 100.0;

        let (mapq_vec, align_score_vec, align_len_vec, as_al_vec, align_prop_vec, align_acc_vec) = 
            if r.is_read1() {
                (&mut self.r1_mapq, &mut self.r1_align_score, &mut self.r1_align_length, &mut self.r1_as_al, &mut self.r1_align_proportion, &mut self.r1_align_accuracy)
            } else {
                (&mut self.r2_mapq, &mut self.r2_align_score, &mut self.r2_align_length, &mut self.r2_as_al, &mut self.r2_align_proportion, &mut self.r2_align_accuracy)
            };

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

    pub fn update_insert_size(&mut self, r1: &BamRecord, r2: &BamRecord, min_mapq: usize, max_ins: i32) {
        if r1.mapq as usize >= min_mapq && r2.mapq as usize >= min_mapq {
            if r1.ref_id == r2.ref_id && r1.ref_id != -1 {
                let (fwd, rev) = if r1.pos <= r2.pos { (r1, r2) } else { (r2, r1) };
                let ref_span: i32 = rev.calculate_ref_span().unwrap_or(0) as i32;
                let insert_size: i32 = (rev.pos + ref_span) - fwd.pos;
                if insert_size > 0 && insert_size <= max_ins {
                    self.pe_insert_size.push(insert_size as f64);
                }
            }
        }
    }
}