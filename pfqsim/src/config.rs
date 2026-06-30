//! n2bio/pfqsim/src/config.rs
//! 

use std::fs::File;
use std::io::{self, Error, ErrorKind};
use std::path::{Path, PathBuf};
use serde::{Deserialize, Serialize};

use crate::cli::AbundanceMode;

// ====================================================================
// Configuration (Inputs)
// ====================================================================

#[derive(Deserialize, Serialize, Debug, Clone)]
pub(crate) struct ConfigRow {
    pub(crate) id: String,
    pub(crate) abundance: f64,
    pub(crate) fasta: PathBuf,
    pub(crate) genome_length: usize,
    pub(crate) model: PathBuf,
    pub(crate) circular: bool,
    pub(crate) sub_rate: f64,
    pub(crate) indel_rate: f64,
    pub(crate) read_length: usize,
    pub(crate) r1_fq: Option<PathBuf>,
    pub(crate) r2_fq: Option<PathBuf>,
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
}

// ====================================================================
// Manifest (Output)
// ====================================================================

#[derive(Serialize, Debug, Clone)]
pub(crate) struct ManifestRow {
    // Fields preserved from the configuration
    pub(crate) id: String,
    pub(crate) abundance: f64,
    pub(crate) fasta: PathBuf,
    pub(crate) genome_length: usize,
    pub(crate) model: PathBuf,
    pub(crate) circular: bool,
    pub(crate) sub_rate: f64,
    pub(crate) indel_rate: f64,
    pub(crate) read_length: usize,
    pub(crate) r1_fq: Option<PathBuf>,
    pub(crate) r2_fq: Option<PathBuf>,
    // The calculated read allocation
    pub(crate) calculated_reads: usize,
}

impl ManifestRow {
    /// Upgrades a raw config row into a manifest plan row by appending the read count
    pub(crate) fn from_config_row(row: ConfigRow, calculated_reads: usize) -> Self {
        Self {
            id: row.id,
            abundance: row.abundance,
            fasta: row.fasta,
            genome_length: row.genome_length,
            model: row.model,
            circular: row.circular,
            sub_rate: row.sub_rate,
            indel_rate: row.indel_rate,
            read_length: row.read_length,
            r1_fq: row.r1_fq,
            r2_fq: row.r2_fq,
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