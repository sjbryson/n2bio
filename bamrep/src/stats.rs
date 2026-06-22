//! n2core/bamrep/src/stats.rs
//! 
use serde::Serialize;
use std::collections::HashMap;

use crate::report::ReportConfig;

// ============================================================================
// Stats
// ============================================================================

// Holds the counts for a single read direction (R1 or R2)
#[derive(Serialize, Default, Clone, Debug)]
pub struct ReadClassStats {
    pub total_reads: usize,
    pub primary_mapped: usize,
    pub primary_mapq: usize,
    pub secondary_mapped: usize,
    pub concordant_mapped: usize,
    pub discordant_mapped: usize,
    pub singletons: usize,
    pub mapq_0: usize,
}

// Holds the combined R1/R2 global stats
#[derive(Serialize, Default, Clone, Debug)]
pub struct GlobalStats {
    pub r1: ReadClassStats,
    pub r2: ReadClassStats,
}

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

        // sample variance
        //let variance_denom = if count > 1.0 { count - 1.0 } else { 1.0 };
        //let variance: f64 = data.iter().map(|&value| {
        //    let diff: f64 = mean - value;
        //    diff * diff
        //}).sum::<f64>() / variance_denom;
        
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