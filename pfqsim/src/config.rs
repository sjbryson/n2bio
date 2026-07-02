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
// Configuration (Inputs)
// ============================================================================

#[derive(Deserialize, Serialize, Debug, Clone)]
pub(crate) struct ConfigRow {
    pub(crate) id: String,
    pub(crate) keyword: String,
    pub(crate) abundance: f64,
    pub(crate) fasta: String,
    pub(crate) model: String,
    pub(crate) sub_rate: f64,
    pub(crate) indel_rate: f64,
    pub(crate) read_length: usize,
    
    #[serde(deserialize_with = "deserialize_flexible_bool")]
    pub(crate) circular: bool,
    
    #[serde(skip_deserializing, default)]
    pub(crate) genome_length: usize,
}

/// Case-insensitive parsing for config boolean declarations 
/// (e.g., 'T', 't', 'TRUE', 'True', 'true', 'F', 'false', etc.)
fn deserialize_flexible_bool<'de, D>(deserializer: D) -> Result<bool, D::Error>
where
    D: Deserializer<'de>,
{
    // Extract the raw field content as a temporary string
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

pub(crate) struct Config {
    pub(crate) rows: Vec<ConfigRow>,
}

impl Config {
    /// Parses a community TSV file into a Config instance
    pub(crate) fn from_tsv<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let file: File = File::open(path)?;
        let mut rdr: csv::Reader<File> = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_reader(file);

        let mut rows: Vec<ConfigRow> = Vec::new();
        for result in rdr.deserialize() {
            let row: ConfigRow = result.map_err(|e| Error::new(ErrorKind::InvalidData, e))?;
            rows.push(row);
        }

        Ok(Self { rows })
    }

    /// Validates all genomes, calculates lengths, and checks circularity constraints.
    pub(crate) fn validate_and_compute_lengths(&mut self) -> io::Result<()> {
        for row in &mut self.rows {
            // Call out to the domain model function to calculate metrics
            let (clean_length, contig_count) = ReferenceGenome::parse_fasta_metrics(&row.fasta)?;

            if clean_length == 0 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Validation Error [{}] -> FASTA file contains 0 valid non-N bases.", row.id),
                ));
            }

            if contig_count > 1 && row.circular {
                println!(
                    "[Validation Warning] Genome '{}' contains multiple contigs ({}) but was marked as circular. \
                     Circularity is only valid for single-contig chromosomes. Forcing circular = false.",
                    row.id, contig_count
                );
                row.circular = false;
            }

            row.genome_length = clean_length;
        }
        
        Ok(())
    }

}

// ============================================================================
// Manifest (Output)
// ============================================================================

#[derive(Serialize, Deserialize, Debug, Clone)]
pub(crate) struct ManifestRow {
    // Fields preserved from the configuration
    pub(crate) id: String,
    pub(crate) keyword: String,
    pub(crate) abundance: f64,
    pub(crate) fasta: String,
    pub(crate) genome_length: usize,
    pub(crate) model: String,
    pub(crate) circular: bool,
    pub(crate) sub_rate: f64,
    pub(crate) indel_rate: f64,
    pub(crate) read_length: usize,
    // The calculated read allocation
    pub(crate) calculated_reads: usize,
}

impl ManifestRow {
    /// Upgrades a raw config row into a manifest plan row by appending the read count
    pub(crate) fn from_config_row(row: ConfigRow, calculated_reads: usize) -> Self {
        Self {
            id: row.id,
            keyword: row.keyword,
            abundance: row.abundance,
            fasta: row.fasta,
            genome_length: row.genome_length,
            model: row.model,
            circular: row.circular,
            sub_rate: row.sub_rate,
            indel_rate: row.indel_rate,
            read_length: row.read_length,
            calculated_reads,
        }
    }
}

pub(crate) struct Manifest {
    pub rows: Vec<ManifestRow>,
}

impl Manifest {
    /// Generates a concrete execution manifest from a configuration map
    pub(crate) fn from_config(config: &Config, total_reads: usize, mode: AbundanceMode) -> Self {
        let mut manifest_rows: Vec<ManifestRow> = Vec::with_capacity(config.rows.len());
        if config.rows.is_empty() {
            return Self { rows: manifest_rows };
        }

        match mode {
            AbundanceMode::ReadFraction => {
                let total_weight: f64 = config.rows.iter().map(|r| r.abundance).sum();
                
                for row in &config.rows {
                    let reads: usize = if total_weight > 0.0 {
                        ((row.abundance / total_weight) * total_reads as f64).round() as usize
                    } else {
                        0
                    };
                    manifest_rows.push(ManifestRow::from_config_row(row.clone(), reads));
                }
            }
            AbundanceMode::CopyFraction => {
                // Read allocation is strictly proportional to (copy_abundance * genome_length)
                let weights: Vec<f64> = config.rows.iter()
                    .map(|r| r.abundance * r.genome_length as f64)
                    .collect();
                let total_weight: f64 = weights.iter().sum();
                
                for (i, row) in config.rows.iter().enumerate() {
                    let reads: usize = if total_weight > 0.0 {
                        ((weights[i] / total_weight) * total_reads as f64).round() as usize
                    } else {
                        0
                    };
                    manifest_rows.push(ManifestRow::from_config_row(row.clone(), reads));
                }
            }
        }

        Self { rows: manifest_rows }
    }

    /// Helper to write the manifest to tsv
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
// Analyze config - based on Manifest - includes "target" field
// ============================================================================

#[derive(Deserialize, Debug, Clone)]
pub(crate) struct AnalyzeRow {
    pub(crate) id: String,      // E.g., "Ecoli_K12" -> matches read header `@{id}:...`
    pub(crate) keyword: String, // Keyword used in generating reads @id:keyword:accession:read_number
    pub(crate) reads: usize,    // The calculated read allocation from compose step
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


