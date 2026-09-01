//! n2bio/src/hist.rs
//! Core histogram structures for empirical distribution accumulation, stats calcs, and sampling.

use serde::{ Serialize, Deserialize };


// ============================================================================
// Histogram
// ============================================================================

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Histogram {
    pub min_val: f64,
    pub max_val: f64,
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
        let total: usize = self.total_count();
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

        let target_rank: f64 = (p / 100.0) * total as f64;
        let mut cumulative_count: f64 = 0.0;

        for (i, &count) in self.counts.iter().enumerate() {
            let prev_cumulative: f64 = cumulative_count;
            cumulative_count += count as f64;

            if cumulative_count >= target_rank {
                // Determine the starting physical value of this bin
                let bin_low_bound: f64 = self.min_val + (i as f64) * self.bin_width;
                if count == 0 { return Some(bin_low_bound); }

                // Linear interpolation: trace exactly where the rank falls inside the bin
                let fraction: f64 = (target_rank - prev_cumulative) / count as f64;
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

    /// Returns the bounding range of observations (ignores empty padding bins).
    /// Returns: Option<(ObservedMin, ObservedMax)>
    pub fn observed_range(&self) -> Option<(f64, f64)> {
        let first_populated_idx: usize = self.counts.iter().position(|&c| c > 0)?;
        let last_populated_idx: usize = self.counts.iter().rposition(|&c| c > 0)?;
        let min_obs: f64 = self.min_val + (first_populated_idx as f64) * self.bin_width;
        let max_obs: f64 = self.min_val + (last_populated_idx as f64 + 1.0) * self.bin_width;
        
        Some((min_obs, max_obs))
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // Helper epsilon for comparing floating-point results
    const EPSILON: f64 = 1e-6;

    #[test]
    fn test_new_histogram_initialization() {
        let hist: Histogram = Histogram::new(0.0, 10.0, 2.0);
        assert_eq!(hist.min_val, 0.0);
        assert_eq!(hist.max_val, 10.0);
        assert_eq!(hist.bin_width, 2.0);
        // ceil((10 - 0) / 2) + 1 = 6 bins
        assert_eq!(hist.counts.len(), 6);
        assert_eq!(hist.total_count(), 0);
        assert!(hist.cdf.is_empty());
    }

    #[test]
    fn test_increment_and_clamping() {
        let mut hist: Histogram = Histogram::new(0.0, 10.0, 2.0);
        
        // Exact midpoints/boundaries
        hist.increment(0.0);  // bin 0
        hist.increment(2.1);  // bin 1
        hist.increment(2.2);  // bin 1
        
        // Out-of-bounds clamping
        hist.increment(-5.0); // Clamped to 0.0 (bin 0)
        hist.increment(15.0); // Clamped to 10.0 (bin 5)

        assert_eq!(hist.counts[0], 2);
        assert_eq!(hist.counts[1], 2);
        assert_eq!(hist.counts[5], 1);
        assert_eq!(hist.total_count(), 5);
    }

    #[test]
    fn test_trim_active_and_empty() {
        let mut hist: Histogram = Histogram::new(0.0, 10.0, 1.0);
        hist.increment(2.0);
        hist.increment(4.0);

        hist.trim();
        // Last active index is at value 4.0 (index 4)
        assert_eq!(hist.counts.len(), 5);
        assert_eq!(hist.max_val, 4.0);

        // Test collapse on completely empty histogram
        let mut empty_hist: Histogram = Histogram::new(0.0, 10.0, 1.0);
        empty_hist.trim();
        assert_eq!(empty_hist.counts.len(), 1);
        assert_eq!(empty_hist.max_val, 0.0);
    }

    #[test]
    fn test_cdf_compilation_and_sampling() {
        let mut hist: Histogram = Histogram::new(0.0, 10.0, 2.0);
        hist.increment(0.0); // bin 0 (weight 1)
        hist.increment(4.0); // bin 2 (weight 3)
        hist.increment(4.0);
        hist.increment(4.0);

        hist.compile_cdf();
        assert_eq!(hist.cdf, vec![0.25, 0.25, 1.0, 1.0, 1.0, 1.0]);

        // Sampling boundary conditions
        assert_eq!(hist.sample(0.0), 0.0);  // First bin
        assert_eq!(hist.sample(0.20), 0.0); // Bin 0
        assert_eq!(hist.sample(0.26), 4.0); // Bin 2
        assert_eq!(hist.sample(1.00), 4.0); // Clamped to highest active bin index
    }

    #[test]
    fn test_mean_mode_observed_range() {
        let mut hist: Histogram = Histogram::new(0.0, 10.0, 2.0);
        
        // bin_index = round((val - min) / width)
        // val = 0.0 -> (0.0/2.0).round() -> bin 0 (midpoint 1.0)
        // val = 4.0 -> (4.0/2.0).round() -> bin 2 (midpoint 5.0)
        hist.increment(0.0); // 1 count in bin 0 (midpoint 1.0)
        hist.increment(4.0); // 3 counts in bin 2 (midpoint 5.0)
        hist.increment(4.0);
        hist.increment(4.0);

        // Mean = (1 * 1.0 + 3 * 5.0) / 4 = 16.0 / 4 = 4.0
        assert_eq!(hist.mean(), Some(4.0));

        // Mode = Midpoint of highest count bin (bin 2 -> 5.0)
        assert_eq!(hist.mode(), Some(5.0));

        // Observed range = bin 0 start (0.0) to bin 2 end (6.0)
        assert_eq!(hist.observed_range(), Some((0.0, 6.0)));
    }

    #[test]
    fn test_percentile_and_median() {
        let mut hist: Histogram = Histogram::new(0.0, 10.0, 1.0);
        
        // val = 0.0 -> (0.0/1.0).round() -> bin 0 (range 0.0 to 1.0)
        for _ in 0..10 {
            hist.increment(0.0);
        }

        // 50th percentile (Median) target_rank = 5.0
        // Interpolation in bin 0: 0.0 + (5.0 / 10.0) * 1.0 = 0.5
        let median = hist.median().unwrap();
        assert!((median - 0.5).abs() < EPSILON);

        // Out of range percentiles return None
        assert_eq!(hist.percentile(-1.0), None);
        assert_eq!(hist.percentile(101.0), None);
    }

    #[test]
    fn test_population_and_sample_stdev() {
        let mut hist: Histogram = Histogram::new(0.0, 10.0, 2.0);
        // Midpoints: bin 0 -> 1.0, bin 2 -> 5.0
        hist.increment(1.0);
        hist.increment(5.0);

        // Mean = 3.0
        // Squared diffs = (1-3)^2 + (5-3)^2 = 4 + 4 = 8
        // Population stdev = sqrt(8 / 2) = sqrt(4) = 2.0
        // Sample stdev     = sqrt(8 / (2 - 1)) = sqrt(8) ≈ 2.828427
        assert!((hist.stdev().unwrap() - 2.0).abs() < EPSILON);
        assert!((hist.sample_stdev().unwrap() - 8.0f64.sqrt()).abs() < EPSILON);

        // Sample stdev returns None for < 2 observations
        let mut single_item_hist: Histogram = Histogram::new(0.0, 10.0, 1.0);
        single_item_hist.increment(1.0);
        assert_eq!(single_item_hist.sample_stdev(), None);
        assert_ne!(single_item_hist.stdev(), None); // Population stdev works for 1 item
    }

    #[test]
    fn test_serde_skips_cdf() {
        let mut hist: Histogram = Histogram::new(0.0, 10.0, 2.0);
        hist.increment(2.0);
        hist.compile_cdf();

        let json: String = serde_json::to_string(&hist).expect("Serialization failed");
        
        // Ensure "cdf" is omitted from JSON payload
        assert!(!json.contains("cdf"));

        let deserialized: Histogram = serde_json::from_str(&json).expect("Deserialization failed");
        assert!(deserialized.cdf.is_empty());
        assert_eq!(deserialized.total_count(), 1);
    }

    #[test]
    fn test_empty_histogram_stat_fallbacks() {
        let hist: Histogram = Histogram::new(0.0, 10.0, 1.0);
        
        assert_eq!(hist.mean(), None);
        assert_eq!(hist.median(), None);
        assert_eq!(hist.mode(), None);
        assert_eq!(hist.stdev(), None);
        assert_eq!(hist.sample_stdev(), None);
        assert_eq!(hist.observed_range(), None);
        assert_eq!(hist.sample(0.5), 0.0); // Defaults to min_val
    }
}