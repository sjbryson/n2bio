//! n2bio/pfqsim/src/config.rs
//! 

#![allow(unused)]

use std::fs::File;
use std::io::{self, Error, ErrorKind};
use std::path::Path;
use serde::{Deserialize, Serialize};
use serde::de::{self, Deserializer};

use crate::cli::AbundanceMode;
use crate::genome::ReferenceGenome;

// ============================================================================
// Compose Config
// ============================================================================

#[derive(Deserialize, Serialize, Debug, Clone)]
pub(crate) struct ComposeTask {
    pub(crate) id: String,
    pub(crate) keyword: String,
    pub(crate) abundance: f64,
    pub(crate) sub_rate: f64,
    pub(crate) indel_rate: f64,
    pub(crate) fasta: String,
    
    #[serde(deserialize_with = "deserialize_flexible_bool")]
    pub(crate) circular: bool,
    
    #[serde(skip_deserializing, default)]
    pub(crate) genome_length: usize,

    #[serde(skip_deserializing, default)]
    pub(crate) calculated_reads: usize,
}

/// Case-insensitive parsing for config boolean declarations 
/// (e.g., 'T', 't', 'TRUE', 'True', 'true', 'F', 'false', etc.)
fn deserialize_flexible_bool<'de, D>(deserializer: D) -> Result<bool, D::Error>
where
    D: Deserializer<'de>,
{
    let s: String = Deserialize::deserialize(deserializer)?;   
    match s.trim().to_uppercase().as_str() {
        "TRUE" | "T" | "YES" | "Y" | "1" => Ok(true),
        "FALSE" | "F" | "NO" | "N" | "0" => Ok(false),
        _ => Err(de::Error::custom(format!(
            "Invalid boolean value '{}'. Expected true/false, T/F, yes/no, or 1/0 (case-insensitive).", 
            s
        ))),
    }
}

pub(crate) struct ComposeConfig {
    pub(crate) rows: Vec<ComposeTask>,
}

impl ComposeConfig {
    /// Parses a community TSV file into a Config instance
    pub(crate) fn from_tsv<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file: File = File::open(path)?;
        let mut rdr: csv::Reader<File> = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_reader(file);

        let mut rows: Vec<ComposeTask> = Vec::new();
        for result in rdr.deserialize() {
            let row: ComposeTask = result.map_err(|e| Error::new(ErrorKind::InvalidData, e))?;
            rows.push(row);
        }

        Ok(Self { rows })
    }

    /// Validates all genomes, calculates lengths, and checks circularity constraints.
    pub(crate) fn validate_and_compute_lengths(&mut self) -> io::Result<()> {
        for task in &mut self.rows {
            // Call out to the domain model function to calculate metrics
            let (clean_length, contig_count) = ReferenceGenome::parse_fasta_metrics(&task.fasta)?;

            if clean_length == 0 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Validation Error [{}] -> FASTA file contains 0 valid non-N bases.", task.id),
                ));
            }

            if contig_count > 1 && task.circular {
                println!(
                    "[Validation Warning] Genome '{}' contains multiple contigs ({}) but was marked as circular. \
                     Circularity is only valid for single-contig chromosomes. Forcing circular = false.",
                    task.id, contig_count
                );
                task.circular = false;
            }

            task.genome_length = clean_length;
        }
        
        Ok(())
    }


    /// Generates a concrete execution manifest from a configuration map
    pub(crate) fn compute_distributions(&mut self, total_reads: usize, mode: AbundanceMode) {
        if self.rows.is_empty() {
            return;
        }

        match mode {
            AbundanceMode::ReadFraction => {
                let total_weight: f64 = self.rows.iter().map(|t| t.abundance).sum();
                
                for task in &mut self.rows {
                    task.calculated_reads = if total_weight > 0.0 {
                        ((task.abundance / total_weight) * total_reads as f64).round() as usize
                    } else {
                        0
                    };
                }
            }
            AbundanceMode::CopyFraction => {
                // Read allocation is strictly proportional to (abundance * genome_length)
                let weights: Vec<f64> = self.rows.iter()
                    .map(|t| t.abundance * t.genome_length as f64)
                    .collect();
                let total_weight: f64 = weights.iter().sum();
                
                for (i, task) in self.rows.iter_mut().enumerate() {
                    task.calculated_reads = if total_weight > 0.0 {
                        ((weights[i] / total_weight) * total_reads as f64).round() as usize
                    } else {
                        0
                    };
                }
            }
        }
    }

     /// Helper to write the config to tsv
    pub(crate) fn save_tsv<P: AsRef<Path>>(&self, path: P) -> io::Result<()> {
        let mut wtr: csv::Writer<File> = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_path(path)?;

        for row in &self.rows {
            wtr.serialize(row).map_err(|e| Error::new(ErrorKind::Other, e))?;
        }
        wtr.flush()?;
        Ok(())
    }
}

// ============================================================================
// Analyze config - reads manifest.tsv output from pfqsim analyze
// ============================================================================

#[derive(Deserialize, Debug, Clone)]
pub(crate) struct AnalyzeRow {
    pub(crate) id: String,      // E.g., "GRCh38" -> matches read header `@{id}:...`
    pub(crate) keyword: String, // Keyword used in generating reads @id:keyword:accession:read_number
    pub(crate) calculated_reads: usize,    // The calculated read allocation from compose step
}

pub(crate) struct AnalyzeConfig {
    pub(crate) rows: Vec<AnalyzeRow>,
}

impl AnalyzeConfig {
    pub(crate) fn from_tsv<P: AsRef<std::path::Path>>(path: P) -> io::Result<Self> {
        let file: File = File::open(path)?;
        let mut rdr: csv::Reader<File> = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_reader(file);

        let mut rows: Vec<AnalyzeRow> = Vec::new();
        for result in rdr.deserialize() {
            let row: AnalyzeRow = result.map_err(|e| Error::new(ErrorKind::InvalidData, e))?;
            rows.push(row);
        }
        Ok(Self { rows })
    }
}

// ============================================================================
// Analyze config - reads manifest.tsv output from pfqsim analyze
// ============================================================================

#[derive(Deserialize, Debug, Clone)]
pub(crate) struct CompareRow {
    pub(crate) id: String,      // Identifier for an analysis report - used as plot labels
    pub(crate) report: String, // Path to the analysis json report associated with the id
}

pub(crate) struct CompareConfig {
    pub(crate) rows: Vec<CompareRow>,
}

impl CompareConfig {
    pub(crate) fn from_tsv<P: AsRef<std::path::Path>>(path: P) -> io::Result<Self> {
        let file: File = File::open(path)?;
        let mut rdr: csv::Reader<File> = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_reader(file);

        let mut rows: Vec<CompareRow> = Vec::new();
        for result in rdr.deserialize() {
            let row: CompareRow = result.map_err(|e| Error::new(ErrorKind::InvalidData, e))?;
            rows.push(row);
        }
        Ok(Self { rows })
    }
}
