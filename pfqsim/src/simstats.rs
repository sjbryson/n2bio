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

/// Accumulates insert size frequencies during BAM/SAM profiling.
#[derive(Debug, Default)]
pub(crate) struct InsertProfiler {
    pub(crate) counts: BTreeMap<i32, usize>,
    pub(crate) total_modeled: usize,
}

impl InsertProfiler {
    
    /// Records a valid insert size observation
    pub(crate) fn update(&mut self, size: i32) {
        *self.counts.entry(size).or_insert(0) += 1;
        self.total_modeled += 1;
    }

    /// Compiles the accumulated frequencies into a finalized, static InsertModel
    pub(crate) fn build(self) -> InsertModel {
        InsertModel {
            insert_dist: NormalDistParams::from_frequency_map(&self.counts),
        }
    }
}

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

#[derive(Debug, Default)]
pub(crate) struct QualityProfiler {
    pub(crate) r1_counts: Vec<BTreeMap<u8, usize>>,
    pub(crate) r2_counts: Vec<BTreeMap<u8, usize>>,
}

impl QualityProfiler {
    
    /// Updates the positional Q-score distribution for a given read slice (1 or 2).
    pub(crate) fn update(&mut self, qual: &[u8], read: u8) {
        let counts: &mut Vec<BTreeMap<u8, usize>> = match read {
            1 => &mut self.r1_counts,
            2 => &mut self.r2_counts,
            _ => panic!("Read argument must be 1 or 2"),
        };

        if counts.len() < qual.len() {
            counts.resize_with(qual.len(), BTreeMap::new);
        }
        
        for (pos, &q) in qual.iter().enumerate() {
            *counts[pos].entry(q).or_insert(0) += 1;
        }
    }

    /// Compiles accumulated frequencies into a QualityModel, padding or 
    /// truncating to match the target read length exactly.
    pub(crate) fn build(self, read_length: usize) -> QualityModel {
        let process_q = |counts: &[BTreeMap<u8, usize>]| -> Vec<NormalDistParams> {
            let mut params: Vec<NormalDistParams> = Vec::with_capacity(read_length);
            for i in 0..read_length {
                if i < counts.len() {
                    params.push(NormalDistParams::from_frequency_map(&counts[i]));
                } else {
                    // Padding fallback if the reference BAM read length was shorter than requested
                    let last_known = params.last().cloned().unwrap_or(NormalDistParams { mean: 30.0, std_dev: 2.0 });
                    params.push(last_known);
                }
            }
            params
        };

        QualityModel {
            r1_quals: process_q(&self.r1_counts),
            r2_quals: process_q(&self.r2_counts),
        }
    }
}

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

#[derive(Debug, Default)]
pub(crate) struct LibraryProfiler {
    pub(crate) quality: QualityProfiler,
    pub(crate) insert_size: InsertProfiler,
}

impl LibraryProfiler {
    pub(crate) fn new() -> Self {
        Self::default()
    }

    /// Records quality scores for a pair of reads
    pub(crate) fn add_qualities(&mut self, r1_qual: &[u8], r2_qual: &[u8]) {
        self.quality.update(r1_qual, 1);
        self.quality.update(r2_qual, 2);
    }

    /// Records a valid insert size observation
    pub(crate) fn add_insert_size(&mut self, size: i32) {
        self.insert_size.update(size);
    }

    /// Finalizes and compiles all accumulated data into a LibraryModel
    pub(crate) fn build(self, read_length: usize) -> LibraryModel {
        LibraryModel {
            insert_size: self.insert_size.build(),
            quality: self.quality.build(read_length),
        }
    }
}

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
