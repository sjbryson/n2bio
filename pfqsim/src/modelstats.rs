//! n2bio/pfqsim/src/hist.rs
//! Payload structures for empirical sequencing data modeling.

use serde::{Serialize, Deserialize};

use n2core::hist::Histogram;

// ============================================================================
// ModelStats
// ============================================================================

#[derive(Serialize, Deserialize, Debug, Clone)]
pub(crate) struct ModelStats {
    pub(crate) max_read_length: usize,
    pub(crate) insert_sizes: Histogram,
    pub(crate) r1_lengths: Histogram,
    pub(crate) r2_lengths: Histogram,
    pub(crate) r1_qualities: Vec<Histogram>,
    pub(crate) r2_qualities: Vec<Histogram>,
}

impl ModelStats {
    /// Initializes a new ModelStats payload with appropriately sized histograms.
    /// `max_insert_size` bounds the max observable insert size (e.g., 2000.0).
    /// `max_read_length` bounds the length histograms and determines the number of quality histograms.
    pub(crate) fn new(max_read_length: usize, max_insert_size: f64) -> Self {
                
        let mut r1_quals: Vec<Histogram> = Vec::with_capacity(max_read_length);
        let mut r2_quals: Vec<Histogram> = Vec::with_capacity(max_read_length);
        
        for _ in 0..max_read_length {
            r1_quals.push(Histogram::new(0.0, 60.0, 1.0));
            r2_quals.push(Histogram::new(0.0, 60.0, 1.0));
        }

        Self {
            max_read_length,
            insert_sizes: Histogram::new(0.0,  max_insert_size, 1.0),
            r1_lengths:   Histogram::new(0.0, max_read_length as f64, 1.0),
            r2_lengths:   Histogram::new(0.0, max_read_length as f64, 1.0),
            r1_qualities: r1_quals,
            r2_qualities: r2_quals,
        }
    }

    /// Cascades the trim operation to all underlying histograms to strip empty right-side bins.
    /// Reduces JSON payload size before saving the model.
    pub(crate) fn trim(&mut self) {
        self.insert_sizes.trim();
        self.r1_lengths.trim();
        self.r2_lengths.trim();
        
        for h in &mut self.r1_qualities {
            h.trim();
        }
        for h in &mut self.r2_qualities {
            h.trim();
        }
    }

    /// Cascades the CDF compilation to all histograms.
    /// Must be called after loading a model JSON to prepare it for random sampling.
    pub(crate) fn compile_cdf(&mut self) {
        self.insert_sizes.compile_cdf();
        self.r1_lengths.compile_cdf();
        self.r2_lengths.compile_cdf();
        
        for h in &mut self.r1_qualities {
            h.compile_cdf();
        }
        for h in &mut self.r2_qualities {
            h.compile_cdf();
        }
    }
}