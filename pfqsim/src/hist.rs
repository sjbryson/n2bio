//! n2bio/pfqsim/src/hist.rs
//! Core histogram structures for empirical distribution accumulation and sampling.

use serde::{Serialize, Deserialize};


// ============================================================================
// Histogram
// ============================================================================

#[derive(Serialize, Deserialize, Debug, Clone)]
pub(crate) struct Histogram {
    pub(crate) min_val: f64,
    pub(crate) max_val: f64,
    pub(crate) bin_width: f64,
    pub(crate) counts: Vec<usize>,
    // Excluded from JSON payload; built in-memory when loaded for sampling
    #[serde(skip)]
    pub(crate) cdf: Vec<f64>,
}

impl Histogram {
    /// Initializes a new histogram with specified boundaries and bin width.
    pub(crate) fn new(min_val: f64, max_val: f64, bin_width: f64) -> Self {
        let num_bins: usize = ((max_val - min_val) / bin_width).ceil() as usize + 1;
        
        Self { 
            min_val, 
            max_val, 
            bin_width, 
            counts: vec![0; num_bins],
            cdf: Vec::new(),
        }
    }

    /// Increments the bin corresponding to the provided value. 
    /// Values outside the min/max bounds are clamped to the nearest edge bin.
    pub(crate) fn increment(&mut self, value: f64) {
        let clamped: f64 = value.clamp(self.min_val, self.max_val);
        let bin_index: usize = ((clamped - self.min_val) / self.bin_width).round() as usize;
        if let Some(bin) = self.counts.get_mut(bin_index) {
            *bin += 1;
        }
    }
    
    /// Returns the total number of items recorded in this histogram.
    pub(crate) fn total_count(&self) -> usize {
        self.counts.iter().sum()
    }

    /// Removes trailing empty bins from the right side of the histogram 
    /// and adjusts the `max_val` to reduce JSON payload sizes.
    pub(crate) fn trim(&mut self) {
        if let Some(last_active_index) = self.counts.iter().rposition(|&count| count > 0) {
            self.counts.truncate(last_active_index + 1);
            self.max_val = self.min_val + (last_active_index as f64 * self.bin_width);
        } else {
            // Collapse completely empty histograms
            self.counts.truncate(1);
            self.max_val = self.min_val;
        }
    }

    /// Prepares the histogram for random sampling by building a Cumulative Distribution Function (CDF).
    /// This should be called once after loading a pre-populated histogram.
    pub(crate) fn compile_cdf(&mut self) {
        let total: f64 = self.total_count() as f64;
        if total == 0.0 {
            self.cdf = Vec::new();
            return;
        }
        
        let mut cumulative: f64 = 0.0;
        self.cdf = Vec::with_capacity(self.counts.len());
        
        for &count in &self.counts {
            cumulative += count as f64;
            self.cdf.push(cumulative / total);
        }
    }

    /// Draws a random value from the empirical distribution. 
    /// `p` should be a uniformly distributed random float in the range [0.0, 1.0].
    pub(crate) fn sample(&self, p: f64) -> f64 {
        if self.cdf.is_empty() {
            return self.min_val; // Safe fallback if sampling an empty distribution
        }
        
        // Find the first index where CDF >= p
        let idx: usize = self.cdf.partition_point(|&prob| prob < p);
        let safe_idx: usize = idx.min(self.counts.len().saturating_sub(1));
        
        self.min_val + (safe_idx as f64 * self.bin_width)
    }

}
