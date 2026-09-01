//! peat/src/bamstats.rs
//! 
use serde::Serialize;
use std::collections::HashMap;

use n2bio::bam::{ BamRecord, BamFlags, BamStats };
use n2bio::hist::Histogram;

// ============================================================================
// Global Stats
// ============================================================================

#[derive(Serialize, Default, Clone, Debug)]
pub struct ReadStats {
    pub total_reads: usize,
    pub total_unaligned: usize,
    pub total_alignments: usize,
    pub primary_mapped: usize,
    pub primary_mapq: usize,
    pub concordant_mapped: usize,
    pub discordant_mapped: usize,
    pub singletons: usize,
    pub secondary_mapped: usize,
    pub mapq_0: usize,
}

#[derive(Serialize, Default, Clone, Debug)]
pub struct GlobalStats {
    pub r1: ReadStats,
    pub r2: ReadStats,
}

impl GlobalStats {
    pub fn update_counts(&mut self, r: &BamRecord) {
        let mate: &mut ReadStats = if r.is_read1() {
            &mut self.r1
        } else if r.is_read2() {
            &mut self.r2
        } else {
            return;
        };
        
        if !r.is_mapped() { mate.total_unaligned += 1; }
        if r.is_mapped() {
            mate.total_alignments += 1;
            
            if r.is_primary() {
                mate.primary_mapped += 1;
                if r.mapq > 0 { mate.primary_mapq += 1; }
                if r.mapq == 0 { mate.mapq_0 += 1;}
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
        self.r1.total_reads = n;
        self.r2.total_reads = n;
    }
}

// ============================================================================
// Summary Stats
// ============================================================================

#[derive(Serialize)]
pub struct ReportData {
    pub global_stats: GlobalStats,
    pub alignment_stats: HashMap<String, StatSummary>,
}

#[derive(Serialize)]
pub struct StatSummary {
    pub count: f64,
    pub mean: f64,
    pub median: f64,
    pub stdev: f64,
    pub min: f64,
    pub max: f64,
    pub histogram: Histogram,
}

impl Default for StatSummary {
    fn default() -> Self {
        Self {
            count: 0.0,
            mean: 0.0,
            median: 0.0,
            stdev: 0.0,
            min: 0.0,
            max: 0.0,
            histogram: Histogram::new(0.0, 0.0, 1.0),
        }
    }
}

impl StatSummary {
    /// Extracts stat summaries from histogram.
    pub fn from_histogram(mut hist: Histogram) -> Self {
        if hist.total_count() == 0 {
            return Self::default();
        }

        // Drop empty trailing bins before sending to JSON
        hist.trim();

        let (min, max) = hist.observed_range().unwrap_or((hist.min_val, hist.max_val));

        StatSummary {
            count: hist.total_count() as f64,
            mean: hist.mean().unwrap_or(0.0),
            median: hist.median().unwrap_or(0.0),
            stdev: hist.stdev().unwrap_or(0.0),
            min,
            max,
            histogram: hist,
        }
    }
}

// ============================================================================
// Stats Accumulator
// ============================================================================

pub struct StatsAccumulator {
    pub pe_insert_size: Histogram,

    pub r1_mapq: Histogram,
    pub r1_align_score: Histogram,
    pub r1_align_length: Histogram,
    pub r1_base_score: Histogram,
    pub r1_align_proportion: Histogram,
    pub r1_align_identity: Histogram,

    pub r2_mapq: Histogram,
    pub r2_align_score: Histogram,
    pub r2_align_length: Histogram,
    pub r2_base_score: Histogram,
    pub r2_align_proportion: Histogram,
    pub r2_align_identity: Histogram,
}

impl StatsAccumulator {
    /// Initialize all histograms
    pub fn new(max_len: usize, max_ins: i32) -> Self {
        let max_al: f64 = max_len as f64;
        let max_as: f64 = 2.0 * max_al;
        let max_ins_f: f64 = max_ins as f64;
      
        Self {
            pe_insert_size: Histogram::new(0.0, max_ins_f, 10.0),

            r1_mapq: Histogram::new(0.0, 60.0, 1.0),
            r1_align_score: Histogram::new(0.0, max_as, 1.0),
            r1_align_length: Histogram::new(0.0, max_al, 1.0),
            r1_base_score: Histogram::new(0.0, 2.0, 0.01),
            r1_align_proportion: Histogram::new(0.0, 1.0, 0.01),
            r1_align_identity: Histogram::new(0.0, 100.0, 0.25),

            r2_mapq: Histogram::new(0.0, 60.0, 1.0),
            r2_align_score: Histogram::new(0.0, max_as, 1.0),
            r2_align_length: Histogram::new(0.0, max_al, 1.0),
            r2_base_score: Histogram::new(0.0, 2.0, 0.01),
            r2_align_proportion: Histogram::new(0.0, 1.0, 0.01),
            r2_align_identity: Histogram::new(0.0, 100.0, 0.25),
        }
    }

    pub fn update_read_stats(&mut self, r: &BamRecord) {
        if r.mapq == 0 || !r.is_primary() || !r.is_mapped() {
            return;
        }

        let (mapq_hist, as_hist, al_hist, bs_hist, ap_hist, ai_hist) = 
            if r.is_read1() {
                (&mut self.r1_mapq, &mut self.r1_align_score, &mut self.r1_align_length, &mut self.r1_base_score, &mut self.r1_align_proportion, &mut self.r1_align_identity)
            } else {
                (&mut self.r2_mapq, &mut self.r2_align_score, &mut self.r2_align_length, &mut self.r2_base_score, &mut self.r2_align_proportion, &mut self.r2_align_identity)
            };

        mapq_hist.increment(r.mapq as f64);
        
        if let Some(val) = r.get_int_tag(b"AS") { 
            as_hist.increment(val as f64); 
        }
        if let Some(val) = r.calculate_alignment_length() { 
            al_hist.increment(val as f64); 
        }
        if let Some(val) = r.calculate_base_score() { 
            bs_hist.increment(val as f64); 
        }
        if let Some(val) = r.calculate_alignment_proportion() { 
            ap_hist.increment(val as f64); 
        }
        if let Some(val) = r.calculate_alignment_identity() { 
            ai_hist.increment(val as f64); 
        }
    }

    pub fn update_insert_size(&mut self, r1: &BamRecord, r2: &BamRecord, min_mapq: usize, max_ins: i32) {
        if r1.mapq as usize >= min_mapq && r2.mapq as usize >= min_mapq {
            if r1.ref_id == r2.ref_id && r1.ref_id != -1 {
                let (fwd, rev) = if r1.pos <= r2.pos { (r1, r2) } else { (r2, r1) };
                let ref_span: i32 = rev.calculate_ref_span().unwrap_or(0) as i32;
                let insert_size: i32 = (rev.pos + ref_span) - fwd.pos;
                
                // > 0 prevents logging overlapping read weirdness/errors as negative insert sizes
                if insert_size > 0 {
                    self.pe_insert_size.increment(insert_size as f64);
                }
            }
        }
    }
}