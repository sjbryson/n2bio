//! n2bio/n2core/src/hist.rs
//! Core histogram structures for empirical distribution accumulation and sampling.

use serde::{Serialize, Deserialize};


// ============================================================================
// Histogram
// ============================================================================

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Histogram {
    pub(crate) min_val: f64,
    pub(crate) max_val: f64,
    pub bin_width: f64,
    pub counts: Vec<usize>,
    // Excluded from JSON payload; built in-memory when loaded for sampling
    #[serde(skip)]
    pub cdf: Vec<f64>,
}

impl Histogram {
    /// Initializes a new histogram with specified boundaries and bin width.
    pub fn new(min_val: f64, max_val: f64, bin_width: f64) -> Self {
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
    pub fn increment(&mut self, value: f64) {
        let clamped: f64 = value.clamp(self.min_val, self.max_val);
        let bin_index: usize = ((clamped - self.min_val) / self.bin_width).round() as usize;
        if let Some(bin) = self.counts.get_mut(bin_index) {
            *bin += 1;
        }
    }
    
    /// Returns the total number of items recorded in this histogram.
    pub fn total_count(&self) -> usize {
        self.counts.iter().sum()
    }

    /// Removes trailing empty bins from the right side of the histogram 
    /// and adjusts the `max_val` to reduce JSON payload sizes.
    pub fn trim(&mut self) {
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
    pub fn compile_cdf(&mut self) {
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
    pub fn sample(&self, p: f64) -> f64 {
        if self.cdf.is_empty() {
            return self.min_val; // Safe fallback if sampling an empty distribution
        }
        
        // Find the first index where CDF >= p
        let idx: usize = self.cdf.partition_point(|&prob| prob < p);
        let safe_idx: usize = idx.min(self.counts.len().saturating_sub(1));
        
        self.min_val + (safe_idx as f64 * self.bin_width)
    }


    /// Helper to get the mathematical center (midpoint) of a specific bin index.
    pub fn bin_midpoint(&self, index: usize) -> f64 {
        self.min_val + (index as f64 + 0.5) * self.bin_width
    }

    /// Calculates the weighted mean of the binned data.
    pub fn mean(&self) -> Option<f64> {
        let total = self.total_count();
        if total == 0 { return None; }

        let weighted_sum: f64 = self.counts.iter().enumerate()
            .map(|(i, &count)| count as f64 * self.bin_midpoint(i))
            .sum();

        Some(weighted_sum / total as f64)
    }

    /// Calculates a precise percentile using linear interpolation across the target bin.
    /// Pass 50.0 for the Median, 25.0 for Q1, 75.0 for Q3, etc.
    pub fn percentile(&self, p: f64) -> Option<f64> {
        if p < 0.0 || p > 100.0 { return None; }
        let total = self.total_count();
        if total == 0 { return None; }

        let target_rank = (p / 100.0) * total as f64;
        let mut cumulative_count = 0.0;

        for (i, &count) in self.counts.iter().enumerate() {
            let prev_cumulative = cumulative_count;
            cumulative_count += count as f64;

            if cumulative_count >= target_rank {
                // Determine the starting physical value of this bin
                let bin_low_bound = self.min_val + (i as f64) * self.bin_width;
                if count == 0 { return Some(bin_low_bound); }

                // Linear interpolation: trace exactly where the rank falls inside the bin
                let fraction = (target_rank - prev_cumulative) / count as f64;
                return Some(bin_low_bound + (fraction * self.bin_width));
            }
        }
        Some(self.max_val)
    }

    /// Convenience wrapper to return the 50th percentile.
    pub fn median(&self) -> Option<f64> {
        self.percentile(50.0)
    }

    /// Calculates the descriptive standard deviation of the specific data points in the histogram
    /// (using N as the denominator to evaluate spread of data).
    pub fn stdev(&self) -> Option<f64> {
        let total: usize = self.total_count();
        if total == 0 { return None; } 
        let mean_val: f64 = self.mean()?;

        let squared_diff_sum: f64 = self.counts.iter().enumerate()
            .map(|(i, &count)| {
                let diff: f64 = self.bin_midpoint(i) - mean_val;
                count as f64 * diff * diff
            })
            .sum();

        Some((squared_diff_sum / total as f64).sqrt())
    }

    /// Calculates the sample standard deviation ($\sigma$).
    pub fn sample_stdev(&self) -> Option<f64> {
        let total: usize = self.total_count();
        // Requires at least 2 observations for a sample standard deviation variance denominator (N - 1)
        if total < 2 { return None; } 
        let mean_val: f64 = self.mean()?;

        let squared_diff_sum: f64 = self.counts.iter().enumerate()
            .map(|(i, &count)| {
                let diff: f64 = self.bin_midpoint(i) - mean_val;
                count as f64 * diff * diff
            })
            .sum();

        Some((squared_diff_sum / (total - 1) as f64).sqrt())
    }

    /// Returns the midpoint of the highest-frequency bin (the peak).
    pub fn mode(&self) -> Option<f64> {
        if self.total_count() == 0 { return None; }
        
        let (max_idx, &max_count) = self.counts.iter().enumerate()
            .max_by_key(|&(_, &count)| count)?;

        if max_count == 0 { return None; }

        Some(self.bin_midpoint(max_idx))
    }

    /// Returns the true bounding range of actual observations (ignores empty padding bins).
    /// Returns: Option<(ObservedMin, ObservedMax)>
    pub fn observed_range(&self) -> Option<(f64, f64)> {
        let first_populated_idx: usize = self.counts.iter().position(|&c| c > 0)?;
        let last_populated_idx: usize = self.counts.iter().rposition(|&c| c > 0)?;
        let min_obs: f64 = self.min_val + (first_populated_idx as f64) * self.bin_width;
        let max_obs: f64 = self.min_val + (last_populated_idx as f64 + 1.0) * self.bin_width;
        
        Some((min_obs, max_obs))
    }
}
