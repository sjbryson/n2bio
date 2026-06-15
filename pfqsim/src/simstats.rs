//! n2bio/pfqsim/src/simstats.rs
//! 

use std::collections::BTreeMap;
use serde::{Deserialize, Serialize};

// ============================================================================
// Core Statistical Parameters
// ============================================================================

/// Represents a standard normal distribution profile
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NormalDistParams {
    pub mean: f64,
    pub std_dev: f64,
}

impl NormalDistParams {
    /// Converts a frequency map into Normal Distribution parameters
    pub fn from_frequency_map<T>(distribution: &BTreeMap<T, usize>) -> Self 
    where
        T: Copy + Into<f64>,
    {
        let total_count: usize = distribution.values().sum();
        
        if total_count == 0 {
            return Self { mean: 0.0, std_dev: 0.0 };
        }

        let n_f64: f64 = total_count as f64;

        // 1. Calculate Mean
        let mut sum: f64 = 0.0;
        for (&value, &count) in distribution {
            sum += value.into() * (count as f64);
        }
        let mean: f64 = sum / n_f64;

        // 2. Calculate Standard Deviation
        let mut variance_sum: f64 = 0.0;
        for (&value, &count) in distribution {
            let diff: f64 = value.into() - mean;
            variance_sum += (count as f64) * diff * diff;
        }
        
        let std_dev = if total_count > 1 {
            (variance_sum / n_f64).sqrt()
        } else {
            0.0
        };

        Self { mean, std_dev }
    }
}

// ============================================================================
// Insert size model
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InsertModel {
    pub insert_dist: NormalDistParams,
}

// ============================================================================
// Quality score model
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityModel {
    pub r1_quals: Vec<NormalDistParams>,
    pub r2_quals: Vec<NormalDistParams>,
}

// ============================================================================
// Library model
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LibraryModel {
    pub insert_size: InsertModel,
    pub quality: QualityModel,
}

impl LibraryModel {
    /// Helper to load a model directly from a JSON file path
    pub fn from_file(path: &str) -> Result<Self, Box<dyn std::error::Error>> {
        let file: std::fs::File = std::fs::File::open(path)?;
        let reader: std::io::BufReader<std::fs::File> = std::io::BufReader::new(file);
        let model: LibraryModel = serde_json::from_reader(reader)?;
        Ok(model)
    }

    /// Serializes and saves the model as a pretty-printed JSON file
    pub fn to_file(&self, path: &str) -> Result<(), Box<dyn std::error::Error>> {
        let file: std::fs::File = std::fs::File::create(path)?;
        let writer: std::io::BufWriter<std::fs::File> = std::io::BufWriter::new(file);
        serde_json::to_writer_pretty(writer, self)?;
        Ok(())
    }
}

// ============================================================================
// Helper functions
// ============================================================================

/// Updates the positional Q-score distribution for a given read.
/// Used during BAM/SAM profiling to accumulate base frequencies.
pub fn update_qscore_model(model: &mut Vec<BTreeMap<u8, usize>>, qual: &[u8]) {
    if model.len() < qual.len() {
        model.resize_with(qual.len(), BTreeMap::new);
    }
    for (pos, &q) in qual.iter().enumerate() {
        *model[pos].entry(q).or_insert(0) += 1;
    }
}