//! n2core/src/fasta.rs
//! FastaRecord struct along with readers and writers.

use std::io::{self, BufRead, Write};
use crate::readers::ReaderType;
use crate::writers::WriterType;

const BUFFER_CAPACITY: usize = 64 * 1024; // 64 KB buffer

// ============================================================================
// FastaRecord struct
// ============================================================================

#[derive(Debug, Clone)]
pub struct FastaRecord {
    pub id:   String,
    pub desc: Option<String>,
    pub seq:  String,
}

impl FastaRecord {
    /// Helper to get the length of the sequence
    pub fn len(&self) -> usize {
        self.seq.len()
    }

    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }
}

// ============================================================================
// FastaReader
// ============================================================================

pub struct FastaReader<R: BufRead> {
    pub reader:  R,
    line_buffer: String,
    next_header: Option<String>,
}

impl<R: BufRead> FastaReader<R> {
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            line_buffer: String::with_capacity(BUFFER_CAPACITY),
            next_header: None,
        }
    }
}

impl FastaReader<ReaderType> {
    /// Automatically detects file type (plain, gz, bz2, stdin) 
    /// based on the path extension.
    pub fn open(path: &str) -> io::Result<Self> {
        Ok(Self::new(ReaderType::open(path)?))
    }

    pub fn from_stdin() -> Self {
        Self::new(ReaderType::from_stdin())
    }

    pub fn from_file(path: &str) -> io::Result<Self> {
        Ok(Self::new(ReaderType::from_file(path)?))
    }

    pub fn from_gz(path: &str) -> io::Result<Self> {
        Ok(Self::new(ReaderType::from_gz(path)?))
    }
}

impl<R: BufRead> Iterator for FastaReader<R> {
    type Item = io::Result<FastaRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        // 1. Determine the header for the current record
        let header_line: String = match self.next_header.take() {
            Some(h) => h,
            None => {
                loop {
                    self.line_buffer.clear();
                    match self.reader.read_line(&mut self.line_buffer) {
                        Ok(0) => return None, // EOF
                        Ok(_) => {
                            let trimmed: &str = self.line_buffer.trim_end();
                            if trimmed.starts_with('>') {
                                break trimmed.to_string();
                            }
                        }
                        Err(e) => return Some(Err(e)),
                    }
                }
            }
        };

        // 2. Parse the ID and Description from the header
        let (id, desc) = parse_fasta_header(&header_line);
        
        // 3. Read the sequence lines until we hit the next '>' or EOF
        let mut seq: String = String::with_capacity(2048);
        
        loop {
            self.line_buffer.clear();
            match self.reader.read_line(&mut self.line_buffer) {
                Ok(0) => break, // EOF, return what we have
                Ok(_) => {
                    let trimmed: &str = self.line_buffer.trim_end();
                    if trimmed.is_empty() {
                        continue;
                    }
                    if trimmed.starts_with('>') {
                        self.next_header = Some(trimmed.to_string());
                        break;
                    }
                    seq.push_str(trimmed);
                }
                Err(e) => return Some(Err(e)),
            }
        }

        Some(Ok(FastaRecord { id, desc, seq }))
    }
}

fn parse_fasta_header(header: &str) -> (String, Option<String>) {
    let text: &str = header.trim_start_matches('>');
    let mut parts = text.splitn(2, |c: char| c.is_ascii_whitespace());
    let id: String = parts.next().unwrap_or("").to_string();
    let desc: Option<String> = parts.next().map(|s| s.trim().to_string());
    
    (id, desc)
}

// ============================================================================
// FastaWriter
// ============================================================================

pub struct FastaWriter<W: Write> {
    pub writer: W,
}

impl<W: Write> FastaWriter<W> {
    pub fn new(writer: W) -> Self {
        Self { writer }
    }

    /// Writes a record using raw bytes.
    pub fn write_record(&mut self, rec: &FastaRecord) -> io::Result<()> {
        self.writer.write_all(b">")?;
        self.writer.write_all(rec.id.as_bytes())?;
        
        if let Some(desc) = &rec.desc {
            self.writer.write_all(b" ")?;
            self.writer.write_all(desc.as_bytes())?;
        }
        
        self.writer.write_all(b"\n")?;
        self.writer.write_all(rec.seq.as_bytes())?;
        self.writer.write_all(b"\n")
    }

    /// Flushes the underlying writer.
    pub fn flush(&mut self) -> io::Result<()> {
        self.writer.flush()
    }

    /// Consumes the FastaWriter, returning the underlying writer.
    pub fn into_inner(self) -> W {
        self.writer
    }
}

impl FastaWriter<WriterType> {
    /// Automatically detects file type (stdout, plain, gz, bz2) 
    /// based on the path extension.
    pub fn create(path: &str) -> io::Result<Self> {
        Ok(Self::new(WriterType::create(path)?))
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fasta_record_len_and_empty() {
        let full_record: FastaRecord = FastaRecord {
            id: "seq1".to_string(),
            desc: None,
            seq: "ACGTACGT".to_string(),
        };
        assert_eq!(full_record.len(), 8);
        assert!(!full_record.is_empty());

        let empty_record: FastaRecord = FastaRecord {
            id: "seq2".to_string(),
            desc: Some("This is an empty sequence".to_string()),
            seq: "".to_string(),
        };
        assert_eq!(empty_record.len(), 0);
        assert!(empty_record.is_empty());
    }

    #[test]
    fn test_parse_fasta_header() {
        // 1. Standard header with ID and description
        let (id, desc) = parse_fasta_header(">chr1 a description");
        assert_eq!(id, "chr1");
        assert_eq!(desc, Some("a description".to_string()));

        // 2. Header with only an ID
        let (id, desc) = parse_fasta_header(">contig_42");
        assert_eq!(id, "contig_42");
        assert_eq!(desc, None);

        // 3. Header with irregular spacing
        let (id, desc) = parse_fasta_header(">read123    some   extra   spaces  ");
        assert_eq!(id, "read123");
        // It trims the end, but preserves internal spaces in the description
        assert_eq!(desc, Some("some   extra   spaces".to_string()));

        // 4. Edge case: Empty header
        let (id, desc) = parse_fasta_header(">");
        assert_eq!(id, "");
        assert_eq!(desc, None);
    }
}
