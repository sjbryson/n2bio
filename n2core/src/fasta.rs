//! n2core/src/fasta.rs

use std::io::{self, BufRead, Write};
use crate::readers::ReaderType;

#[derive(Debug, Clone)]
pub struct FastaRecord {
    pub id: String,
    pub desc: Option<String>,
    pub seq: String,
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

/// READER
/// 
/// 
/// 
/// 
pub struct FastaReader<R: BufRead> {
    pub reader: R,
    line_buffer: String,
    next_header: Option<String>,
}

impl<R: BufRead> FastaReader<R> {
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            line_buffer: String::with_capacity(1024),
            next_header: None,
        }
    }
}

impl FastaReader<ReaderType> {
    pub fn from_file(path: &str) -> io::Result<Self> {
        Ok(Self::new(ReaderType::from_file(path)?))
    }

    pub fn from_gz(path: &str) -> io::Result<Self> {
        Ok(Self::new(ReaderType::from_gz(path)?))
    }

    pub fn from_stdin() -> Self {
        Self::new(ReaderType::from_stdin())
    }
}

impl<R: BufRead> Iterator for FastaReader<R> {
    type Item = io::Result<FastaRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        // 1. Determine the header for the current record
        let header_line = match self.next_header.take() {
            Some(h) => h,
            None => {
                loop {
                    self.line_buffer.clear();
                    match self.reader.read_line(&mut self.line_buffer) {
                        Ok(0) => return None, // EOF
                        Ok(_) => {
                            let trimmed = self.line_buffer.trim_end();
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
        let mut seq = String::with_capacity(2048);
        
        loop {
            self.line_buffer.clear();
            match self.reader.read_line(&mut self.line_buffer) {
                Ok(0) => break, // EOF, return what we have
                Ok(_) => {
                    let trimmed = self.line_buffer.trim_end();
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
    let text = header.trim_start_matches('>');
    let mut parts = text.splitn(2, |c: char| c.is_ascii_whitespace());
    
    let id = parts.next().unwrap_or("").to_string();
    let desc = parts.next().map(|s| s.trim().to_string());
    
    (id, desc)
}

/// WRITER
/// 
/// 
/// 
/// 
/// 
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