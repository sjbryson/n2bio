//! n2bio/peat/src/binpooler.rs
//! 

use serde::Serialize;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{self, BufWriter};
use std::path::Path;

#[derive(Debug, Default, Serialize)]
pub struct BinStats {
    pub concordant_pairs: usize,
    pub discordant_pairs: usize,
    pub r1_orphans: usize,
    pub r2_orphans: usize,
    pub total_pairs_written: usize,
}

#[derive(Debug, Default, Serialize)]
pub struct BinReadReport {
    pub total_bins: usize,
    pub total_pairs_binned: usize,
    pub bins: BTreeMap<String, BinStats>,
}

impl BinReadReport {
    pub fn inc_concordant(&mut self, bin_id: &str) {
        let entry: &mut BinStats = self.bins.entry(bin_id.to_string()).or_default();
        entry.concordant_pairs += 1;
        entry.total_pairs_written += 1;
        self.total_pairs_binned += 1;
    }

    pub fn inc_discordant(&mut self, bin_id: &str) {
        let entry: &mut BinStats = self.bins.entry(bin_id.to_string()).or_default();
        entry.discordant_pairs += 1;
        entry.total_pairs_written += 1;
        self.total_pairs_binned += 1;
    }

    pub fn inc_r1_orphan(&mut self, bin_id: &str) {
        let entry: &mut BinStats = self.bins.entry(bin_id.to_string()).or_default();
        entry.r1_orphans += 1;
        entry.total_pairs_written += 1;
        self.total_pairs_binned += 1;
    }

    pub fn inc_r2_orphan(&mut self, bin_id: &str) {
        let entry: &mut BinStats = self.bins.entry(bin_id.to_string()).or_default();
        entry.r2_orphans += 1;
        entry.total_pairs_written += 1;
        self.total_pairs_binned += 1;
    }

    /// Writes the stats report as a pretty-printed JSON file.
    pub fn write_json<P: AsRef<Path>>(&mut self, path: P) -> io::Result<()> {
        self.total_bins = self.bins.len();
        let file: File = File::create(path)?;
        let writer: BufWriter<File> = BufWriter::new(file);
        serde_json::to_writer_pretty(writer, self)?;
        Ok(())
    }
}