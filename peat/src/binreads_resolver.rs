//! n2bio/peat/src/binresolver.rs
//! 

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

// ============================================================================
// BinResolver
// ============================================================================
pub struct BinResolver {
    target_to_bin: Option<HashMap<String, String>>,
}

impl BinResolver {
    /// Constructs a `BinResolver`.
    /// 
    /// If `reference_map` is provided and non-empty, reads the TSV mapping file (target -> bin).
    /// If `reference_map` is empty or `None`, falls back to using reference names directly as bins.
    pub fn new(reference_map: Option<&str>) -> io::Result<Self> {
        if let Some(map_path) = reference_map {
            if !map_path.trim().is_empty() {
                let target_to_bin: HashMap<String, String> = Self::parse_tsv(map_path)?;
                return Ok(Self {
                    target_to_bin: Some(target_to_bin),
                });
            }
        }

        Ok(Self {
            target_to_bin: None,
        })
    }

    /// Helper to parse TSV file: reference_id (col 0) -> bin_id (col 1)
    fn parse_tsv(path_str: &str) -> io::Result<HashMap<String, String>> {
        let path: &Path = Path::new(path_str);
        let file: File = File::open(path)?;
        let reader: BufReader<File> = BufReader::new(file);

        let mut map: HashMap<String, String> = HashMap::new();

        for (line_num, line_result) in reader.lines().enumerate() {
            let line: String = line_result?;
            let trimmed: &str = line.trim();

            // Skip comment lines and empty lines
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }

            let mut fields: std::str::Split<'_, char> = trimmed.split('\t');
            let target_id: Option<&str> = fields.next();
            let bin_id: Option<&str> = fields.next();

            match (target_id, bin_id) {
                (Some(target), Some(bin)) if !target.is_empty() && !bin.is_empty() => {
                    map.insert(target.to_string(), bin.to_string());
                }
                _ => {
                    // Fail early on malformed TSV lines
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "Malformed TSV line {} in {}: expected 2 tab-separated columns (target_id\\tbin_id)",
                            line_num + 1,
                            path_str
                        ),
                    ));
                }
            }
        }

        Ok(map)
    }

    /// Resolves a target reference name to its corresponding bin ID.
    /// 
    /// If a custom mapping was loaded, returns the mapped bin ID (or `None` if unmapped in TSV).
    /// If no mapping was provided, returns the `target_name` itself as the bin ID.
    pub fn get_bin<'a>(&'a self, target_name: &'a str) -> Option<&'a str> {
        match &self.target_to_bin {
            Some(map) => map.get(target_name).map(|s| s.as_str()),
            None => Some(target_name),
        }
    }
}