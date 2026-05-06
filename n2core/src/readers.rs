//! n2core/src/readers.rs
//! Buffered readers for plain text and compressed files.

use std::io::{self, BufReader, Read, BufRead, Stdin};
use std::fs::File;
use flate2::read::MultiGzDecoder;
use bzip2::read::BzDecoder;

const BUFFER_CAPACITY: usize = 64 * 1024; // 64 KB buffer

// ============================================================================
// Enums
// ============================================================================

/// Reader capable of handling standard input, 
/// plain text files, and compressed files (gzip, bzip2).
pub enum ReaderType {
    Stdin(BufReader<Stdin>),
    File(BufReader<File>),
    Gz(BufReader<MultiGzDecoder<File>>),
    Bz2(BufReader<BzDecoder<File>>),
}

// ============================================================================
// Constructors
// ============================================================================

impl ReaderType {
    /// Automatically detects the correct reader type based on the 
    /// file extension. Uses `from_stdin` if the path is `-`.
    pub fn open(path: &str) -> io::Result<Self> {
        if path == "-" {
            Ok(Self::from_stdin())
        } else if path.ends_with(".gz") {
            Self::from_gz(path)
        } else if path.ends_with(".bz2") {
            Self::from_bz(path)
        } else {
            Self::from_file(path)
        }
    }

    /// Read from standard input.
    pub fn from_stdin() -> Self {
        let buffer: BufReader<Stdin> = BufReader::with_capacity(BUFFER_CAPACITY, io::stdin());
        
        ReaderType::Stdin(buffer)
    }

    /// Buffered reader for uncompressed plain text files.
    pub fn from_file(path: &str) -> io::Result<Self> {
        let file: File = File::open(path)?;
        let buffer: BufReader<File> = BufReader::with_capacity(BUFFER_CAPACITY, file);
        
        Ok(ReaderType::File(buffer))
    }

    /// Buffered, gzipped reader capable of handling concatenated gzip files.
    pub fn from_gz(path: &str) -> io::Result<Self> {
        let file: File = File::open(path)?;
        let decoder: MultiGzDecoder<File> = MultiGzDecoder::new(file);
        let buffer: BufReader<MultiGzDecoder<File>> = BufReader::with_capacity(BUFFER_CAPACITY, decoder); 
        
        Ok(ReaderType::Gz(buffer))
    }

    /// Buffered, parallel bzip2 reader.
    pub fn from_bz(path: &str) -> io::Result<Self> {
        let file: File = File::open(path)?;
        let decoder: BzDecoder<File> = BzDecoder::new(file);
        let buffer: BufReader<BzDecoder<File>> = BufReader::with_capacity(BUFFER_CAPACITY, decoder);
        
        Ok(ReaderType::Bz2(buffer))
    }
}

// ============================================================================
// Trait Implementations
// ============================================================================

impl Read for ReaderType {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        match self {
            ReaderType::Stdin(r) => r.read(buf),
            ReaderType::File(r) => r.read(buf),
            ReaderType::Gz(r) => r.read(buf),
            ReaderType::Bz2(r) => r.read(buf),
        }
    }
}

impl BufRead for ReaderType {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        match self {
            ReaderType::Stdin(r) => r.fill_buf(),
            ReaderType::File(r) => r.fill_buf(),
            ReaderType::Gz(r) => r.fill_buf(),
            ReaderType::Bz2(r) => r.fill_buf(),
        }
    }

    fn consume(&mut self, amt: usize) {
        match self {
            ReaderType::Stdin(r) => r.consume(amt),
            ReaderType::File(r) => r.consume(amt),
            ReaderType::Gz(r) => r.consume(amt),
            ReaderType::Bz2(r) => r.consume(amt),
        }
    }

    fn read_line(&mut self, buf: &mut String) -> io::Result<usize> {
        match self {
            ReaderType::Stdin(r) => r.read_line(buf),
            ReaderType::File(r) => r.read_line(buf),
            ReaderType::Gz(r) => r.read_line(buf),
            ReaderType::Bz2(r) => r.read_line(buf),
        }
    }
}