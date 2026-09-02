//! n2bio/src/metadata.rs
//! Core structs and traits for incorporating metadata into reports

use csv::ReaderBuilder;
use serde_json::{Map, Number, Value};
use std::collections::HashMap;
use std::path::Path;
use std::fs::File;
use std::io::{BufRead, BufReader};
use thiserror::Error;

// ============================================================================
// Metadata - Errors
// ============================================================================

#[derive(Error, Debug)]
pub enum MetadataError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("CSV parsing error: {0}")]
    Csv(#[from] csv::Error),

    #[error("JSON parsing error: {0}")]
    Json(#[from] serde_json::Error),

    #[error("Key column '{0}' not found in the dataset")]
    MissingKeyColumn(String),

    #[error("A row is missing a value for the key column")]
    MissingKeyValue,

    #[error("Expected a JSON array at the root of the file")]
    ExpectedJsonArray,
}

// ============================================================================
// Metadata
// ============================================================================


pub struct Metadata {
    store: HashMap<String, Value>,
}

impl Metadata {
    pub fn from_tsv<P: AsRef<Path>>(path: P, key_column: &str) -> Result<Self, MetadataError> {
        let mut rdr: csv::Reader<File> = ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(path)?;

        let headers: csv::StringRecord = rdr.headers()?.clone();
        
        let key_index: usize = headers
            .iter()
            .position(|h| h == key_column)
            .ok_or_else(|| MetadataError::MissingKeyColumn(key_column.to_string()))?;

        let mut store: HashMap<String, Value> = HashMap::new();

        for result in rdr.records() {
            let record: csv::StringRecord = result?;
            
            let key_value: String = record.get(key_index)
                .ok_or(MetadataError::MissingKeyValue)?
                .trim()
                .to_string();

            let mut row_map: Map<String, Value> = Map::new();

            for (header, field) in headers.iter().zip(record.iter()) {
                row_map.insert(
                    header.trim().to_string(), 
                    infer_field_type(field)
                );
            }

            store.insert(key_value, Value::Object(row_map));
        }

        Ok(Metadata { store })
    }

    pub fn from_json<P: AsRef<Path>>(path: P, key_column: &str) -> Result<Self, MetadataError> {
        let file: File = File::open(path)?;
        let reader: BufReader<File> = BufReader::new(file);
        
        let data: Value = serde_json::from_reader(reader)?;
        
        let mut store: HashMap<String, Value> = HashMap::new();

        if let Some(array) = data.as_array() {
            for item in array {
                if let Some(obj) = item.as_object() {
                    if let Some(key_value) = obj.get(key_column) {
                        let key_str: String = match key_value {
                            Value::String(s) => s.clone(),
                            other => other.to_string(),
                        };
                        
                        store.insert(key_str, item.clone());
                    }
                }
            }
        } else {
            return Err(MetadataError::ExpectedJsonArray);
        }

        Ok(Metadata { store })
    }

    pub fn from_jsonl<P: AsRef<Path>>(path: P, key_column: &str) -> Result<Self, MetadataError> {
        let file: File = File::open(path)?;
        let reader: BufReader<File> = BufReader::new(file); 
        let mut store: HashMap<String, Value> = HashMap::new();

        for line_result in reader.lines() {
            let line: String = line_result?;
            
            if line.trim().is_empty() {
                continue;
            }

            let item: Value = serde_json::from_str(&line)?;

            if let Some(obj) = item.as_object() {
                if let Some(key_value) = obj.get(key_column) {
                    let key_str: String = match key_value {
                        Value::String(s) => s.clone(),
                        other => other.to_string(),
                    };
                    store.insert(key_str, item);
                }
            }
        }

        Ok(Metadata { store })
    }
    
    pub fn lookup(&self, key: &str) -> Option<&Value> {
        self.store.get(key)
    }
}

pub fn infer_field_type(field: &str) -> Value {
    let trimmed: &str = field.trim();

    if trimmed.is_empty() || trimmed.eq_ignore_ascii_case("null") {
        return Value::Null;
    }
    
    // case-insensitive booleans
    if trimmed.eq_ignore_ascii_case("true") {
        return Value::Bool(true);
    }
    if trimmed.eq_ignore_ascii_case("false") {
        return Value::Bool(false);
    }

    // integers
    if let Ok(i) = trimmed.parse::<i64>() {
        return Value::Number(i.into());
    }

    // floats
    if let Ok(f) = trimmed.parse::<f64>() {
        if let Some(num) = Number::from_f64(f) {
            return Value::Number(num);
        }
    }

    // strings
    Value::String(field.to_string())
}

// ============================================================================
// Metadata - Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_infer_field_type() {
        assert_eq!(infer_field_type("true"), Value::Bool(true));
        assert_eq!(infer_field_type("False"), Value::Bool(false));
        assert_eq!(infer_field_type("42"), Value::Number(42.into()));
        assert_eq!(infer_field_type("3.14"), Value::Number(Number::from_f64(3.14).unwrap()));
        assert_eq!(infer_field_type("null"), Value::Null);
        assert_eq!(infer_field_type(""), Value::Null);
        assert_eq!(infer_field_type("hello"), Value::String("hello".to_string()));
    }

    #[test]
    fn test_from_tsv() {
        let mut file: NamedTempFile = NamedTempFile::new().unwrap();
        // Create TSV
        writeln!(file, "id\tnum_value\tbool_value\tnotes").unwrap();
        writeln!(file, "val_1\t45\ttrue\t").unwrap(); // Empty notes should be Null
        writeln!(file, "val_2\t3.14\tfalse\tnull").unwrap(); // "null" notes should be Null

        let metadata: Metadata = Metadata::from_tsv(file.path(), "id").unwrap();

        // Check id1
        let id1: &Value = metadata.lookup("val_1").expect("val_1 should exist");
        assert_eq!(id1["num_value"], Value::Number(45.into()));
        assert_eq!(id1["bool_value"], Value::Bool(true));
        assert_eq!(id1["notes"], Value::Null);

        // Check id2
        let id2: &Value = metadata.lookup("val_2").expect("val_2 should exist");
        assert_eq!(id2["num_value"], Value::Number(Number::from_f64(3.14).unwrap()));
        assert_eq!(id2["bool_value"], Value::Bool(false));
    }

    #[test]
    fn test_missing_key_column_error() {
        let mut file: NamedTempFile = NamedTempFile::new().unwrap();
        writeln!(file, "name\tvalue").unwrap();
        writeln!(file, "test\t100").unwrap();

        let result: Result<Metadata, MetadataError> = Metadata::from_tsv(file.path(), "id");
        
        // test returns correct error
        assert!(matches!(result, Err(MetadataError::MissingKeyColumn(_))));
    }

    #[test]
    fn test_from_json_array() {
        let mut file: NamedTempFile = NamedTempFile::new().unwrap();
        let json_data: &str = r#"
        [
            {"id": "val_1", "property": "AA"},
            {"id": "val_2", "property": "BB"}
        ]
        "#;
        write!(file, "{}", json_data).unwrap();

        let metadata: Metadata = Metadata::from_json(file.path(), "id").unwrap();
        
        let id2: &Value = metadata.lookup("val_2").unwrap();
        assert_eq!(id2["property"], Value::String("BB".to_string()));
    }

    #[test]
    fn test_from_json_array_fails_on_object() {
        let mut file: NamedTempFile = NamedTempFile::new().unwrap();
        // Write a dict instead of an array
        let json_data: &str = r#"{"id": "val_1", "property": "AA"}"#;
        write!(file, "{}", json_data).unwrap();

        let result: Result<Metadata, MetadataError> = Metadata::from_json(file.path(), "id");
        
        assert!(matches!(result, Err(MetadataError::ExpectedJsonArray)));
    }

    #[test]
    fn test_from_jsonl() {
        let mut file: NamedTempFile = NamedTempFile::new().unwrap();
        let jsonl_data: &str = "\
            {\"id\": \"val_1\", \"property\": \"AA\"}\n\
            {\"id\": \"val_2\", \"property\": \"BB\"}\n\
        ";
        write!(file, "{}", jsonl_data).unwrap();

        let metadata: Metadata = Metadata::from_jsonl(file.path(), "id").unwrap();
        
        let id1: &Value = metadata.lookup("val_1").unwrap();
        assert_eq!(id1["property"], Value::String("AA".to_string()));
    }
}