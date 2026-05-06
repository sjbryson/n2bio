//! n2core/src/writers.rs
//! Buffered writers for plain text and compressed files.

use std::fs::File;
use std::io::{self, BufWriter, Write, Stdout};
use flate2::write::GzEncoder;
use flate2::Compression as GzCompression;
use bzip2::write::BzEncoder;
use bzip2::Compression as BzCompression;
use gzp::{deflate::Gzip, ZBuilder, ZWriter};

const BUFFER_CAPACITY: usize = 64 * 1024; // 64 KB buffer

// ============================================================================
// Enums
// ============================================================================

/// Writer capable of handling standard output, 
/// plain text files, and compressed files (gzip, bzip2).
pub enum WriterType {
    File(BufWriter<File>),
    Stdout(BufWriter<Stdout>),
    Gz(BufWriter<GzEncoder<File>>),
    Bz(BufWriter<BzEncoder<File>>),
    MultiGz(Box<dyn Write>),
}

// ============================================================================
// Constructors
// ============================================================================


impl WriterType {
    /// Automatically detects the correct writer type based on the 
    /// file extension. Uses `to_stdout` if the path is `-`. Defaults to single-threaded compression.
    pub fn create(path: &str) -> io::Result<Self> {
        if path == "-" {
            Ok(Self::to_stdout())
        } else if path.ends_with(".gz") {
            Self::to_gz(path)
        } else if path.ends_with(".bz2") {
            Self::to_bz(path)
        } else {
            Self::to_file(path)
        }
    }

    /// Wraps standard output in a buffered writer.
    pub fn to_stdout() -> Self {
        let buffer: BufWriter<Stdout> = BufWriter::with_capacity(BUFFER_CAPACITY, io::stdout());
        
        WriterType::Stdout(buffer)
    }

    /// Standard buffered writer for uncompressed plain text.
    pub fn to_file(path: &str) -> io::Result<Self> {
        let file: File = File::create(path)?;
        let buffer: BufWriter<File> = BufWriter::with_capacity(BUFFER_CAPACITY, file);
        
        Ok(WriterType::File(buffer))
    }

    /// Buffered, gzipped writer (single-threaded).
    pub fn to_gz(path: &str) -> io::Result<Self> {
        let file: File = File::create(path)?;
        let gz: GzEncoder<File> = GzEncoder::new(file, GzCompression::default());
        let buffer: BufWriter<GzEncoder<File>> = BufWriter::with_capacity(BUFFER_CAPACITY, gz); 
        
        Ok(WriterType::Gz(buffer))
    }

    /// Buffered, bzip2 writer.
    pub fn to_bz(path: &str) -> io::Result<Self> {
        let file: File = File::create(path)?;
        let bz: BzEncoder<File> = BzEncoder::new(file, BzCompression::default());
        let buffer: BufWriter<BzEncoder<File>> = BufWriter::with_capacity(BUFFER_CAPACITY, bz);
        
        Ok(WriterType::Bz(buffer))
    }

    /// Multi-threaded gzip writer.
    pub fn to_multithreaded_gz(path: &str, threads: usize) -> io::Result<Self> {
        let file: File = File::create(path)?;
        let par_writer: Box<dyn ZWriter> = ZBuilder::<Gzip, _>::new()
            .num_threads(threads)
            .from_writer(file);
            
        Ok(WriterType::MultiGz(Box::new(par_writer)))
    }
}

// ============================================================================
// Trait Implementations
// ============================================================================

impl Write for WriterType {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            WriterType::Stdout(w) => w.write(buf),
            WriterType::File(w) => w.write(buf),
            WriterType::Gz(w) => w.write(buf),
            WriterType::Bz(w) => w.write(buf),
            WriterType::MultiGz(w) => w.write(buf),
        }
    }
    fn flush(&mut self) -> io::Result<()> {
        match self {
            WriterType::Stdout(w) => w.flush(),
            WriterType::File(w) => w.flush(),
            WriterType::Gz(w) => w.flush(),
            WriterType::Bz(w) => w.flush(),
            WriterType::MultiGz(w) => w.flush(),
        }
    }
}
