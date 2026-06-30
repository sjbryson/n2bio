//! n2bio/pfqsim/src/simstats.rs
//! 

use rand::Rng;
use rand_distr::{Normal, Distribution};
use std::collections::BTreeMap;
use serde::{Deserialize, Serialize};

// ============================================================================
// Core Statistical Parameters
// ============================================================================

/// Represents a standard normal distribution profile
#[derive(Debug, Clone, Serialize, Deserialize)]
pub(crate) struct NormalDistParams {
    pub(crate) mean: f64,
    pub(crate) std_dev: f64,
}

impl NormalDistParams {
    /// Converts a frequency map into Normal Distribution parameters
    pub(crate) fn from_frequency_map<T>(distribution: &BTreeMap<T, usize>) -> Self 
    where
        T: Copy + Into<f64>,
    {
        let total_count: usize = distribution.values().sum();
        
        if total_count == 0 {
            return Self { mean: 0.0, std_dev: 0.0 };
        }

        let n_f64: f64 = total_count as f64;

        // Calculate Mean
        let mut sum: f64 = 0.0;
        for (&value, &count) in distribution {
            sum += value.into() * (count as f64);
        }
        let mean: f64 = sum / n_f64;

        // Calculate Standard Deviation
        let mut variance_sum: f64 = 0.0;
        for (&value, &count) in distribution {
            let diff: f64 = value.into() - mean;
            variance_sum += (count as f64) * diff * diff;
        }
        
        let std_dev: f64 = if total_count > 1 {
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
pub(crate) struct InsertModel {
    pub(crate) insert_dist: NormalDistParams,
}

pub(crate) struct InsertSize {
    dist: Normal<f64>,
}

impl InsertSize {
    pub(crate) fn new(params: &InsertModel) -> Result<Self, &'static str> {
        // Access the nested NormalDistParams fields within InsertModel
        let mean: f64 = params.insert_dist.mean;
        let std_dev: f64 = if params.insert_dist.std_dev <= 0.0 { 0.1 } else { params.insert_dist.std_dev };
        
        let dist: Normal<f64> = Normal::new(mean, std_dev)
            .map_err(|_| "Failed to create normal distribution for insert size")?;
        
        Ok(Self { dist })
    }

    /// Samples a random insert size, rounding to the nearest integer
    pub(crate) fn sample<R: Rng>(&self, rng: &mut R) -> usize {
        self.dist.sample(rng).round().max(0.0) as usize
    }
}

// ============================================================================
// Quality score model
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub(crate) struct QualityModel {
    pub(crate) r1_quals: Vec<NormalDistParams>,
    pub(crate) r2_quals: Vec<NormalDistParams>,
}

pub(crate) struct QualityScores {
    pub(crate) r1_qual: Vec<Normal<f64>>,
    pub(crate) r2_qual: Vec<Normal<f64>>,
}

impl QualityScores {
    pub(crate) fn new(model: &QualityModel) -> Result<Self, &'static str> {
        let to_normal = |params: &NormalDistParams| {
            let std_dev: f64 = if params.std_dev <= 0.0 { 0.1 } else { params.std_dev };
            Normal::new(params.mean, std_dev).unwrap()
        };

        let r1_qual: Vec<Normal<f64>> = model.r1_quals.iter().map(to_normal).collect();
        let r2_qual: Vec<Normal<f64>> = model.r2_quals.iter().map(to_normal).collect();

        Ok(Self { r1_qual, r2_qual })
    }

    /// Generates a quality string of ASCII characters for the specified read (1 or 2)
    pub(crate) fn generate<R: Rng>(&self, rng: &mut R, length: usize, read: u8) -> Vec<u8> {
        let distributions: &Vec<Normal<f64>> = match read {
            1 => &self.r1_qual,
            2 => &self.r2_qual,
            _ => panic!("Read argument must be 1 or 2"),
        };

        distributions.iter().take(length).map(|dist| {
            let q: f64 = dist.sample(rng).round();
            let clamped_q: u8 = q.clamp(0.0, 41.0) as u8;
            clamped_q + 33
        }).collect()
    }
}

// ============================================================================
// Library model
// ============================================================================

#[derive(Debug, Clone, Serialize, Deserialize)]
pub(crate) struct LibraryModel {
    pub(crate) insert_size: InsertModel,
    pub(crate) quality: QualityModel,
}

impl LibraryModel {
    /// Helper to load a model directly from a JSON file path
    pub(crate) fn from_file(path: &str) -> Result<Self, Box<dyn std::error::Error>> {
        let file: std::fs::File = std::fs::File::open(path)?;
        let reader: std::io::BufReader<std::fs::File> = std::io::BufReader::new(file);
        let model: LibraryModel = serde_json::from_reader(reader)?;
        Ok(model)
    }

    /// Serializes and saves the model as a pretty-printed JSON file
    pub(crate) fn to_file(&self, path: &str) -> Result<(), Box<dyn std::error::Error>> {
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
pub(crate) fn update_qscore_model(model: &mut Vec<BTreeMap<u8, usize>>, qual: &[u8]) {
    if model.len() < qual.len() {
        model.resize_with(qual.len(), BTreeMap::new);
    }
    for (pos, &q) in qual.iter().enumerate() {
        *model[pos].entry(q).or_insert(0) += 1;
    }
}