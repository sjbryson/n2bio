//! n2core/src/writers.rs

use std::fs::File;
use std::io::{self, BufWriter, Write};
use flate2::write::GzEncoder;
use flate2::Compression as GzCompression;
use bzip2::write::BzEncoder;
use bzip2::Compression as BzCompression;
use gzp::{deflate::Gzip, ZBuilder};

pub enum WriterType {
    File(BufWriter<File>),
    Stdout(io::Stdout),
    Gz(BufWriter<GzEncoder<File>>),
    Bz(BufWriter<BzEncoder<File>>),
    MultiGz(Box<dyn Write>),
}

impl Write for WriterType {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            WriterType::File(w) => w.write(buf),
            WriterType::Stdout(w) => w.write(buf),
            WriterType::Gz(w) => w.write(buf),
            WriterType::Bz(w) => w.write(buf),
            WriterType::MultiGz(w) => w.write(buf),
        }
    }
    fn flush(&mut self) -> io::Result<()> {
        match self {
            WriterType::File(w) => w.flush(),
            WriterType::Stdout(w) => w.flush(),
            WriterType::Gz(w) => w.flush(),
            WriterType::Bz(w) => w.flush(),
            WriterType::MultiGz(w) => w.flush(),
        }
    }
}

impl WriterType {
    
    /// Wraps standard output.
    pub fn to_stdout() -> Self {
        WriterType::Stdout(io::stdout())
    }

    /// Standard buffered writer for uncompressed plain text.
    pub fn to_file(path: &str) -> io::Result<Self> {
        let file = File::create(path)?;
        let buffer = BufWriter::with_capacity(64 * 1024, file);
        
        Ok(WriterType::File(buffer))
    }

    /// Buffered, gzipped writer.
    pub fn to_gz(path: &str) -> io::Result<Self> {
        let file = File::create(path)?;
        let gz = GzEncoder::new(file, GzCompression::default());
        // Buffering BEFORE the encoder is highly recommended for performance
        let buffer = BufWriter::with_capacity(64 * 1024, gz); 
        
        Ok(WriterType::Gz(buffer))
    }

    /// Buffered, bzip2 writer.
    pub fn to_bz(path: &str) -> io::Result<Self> {
        let file = File::create(path)?;
        let bz = BzEncoder::new(file, BzCompression::default());
        let buffer = BufWriter::with_capacity(64 * 1024, bz);
        
        Ok(WriterType::Bz(buffer))
    }

    /// Multi-threaded gzip writer.
    pub fn to_multithreaded_gz(path: &str, threads: usize) -> io::Result<Self> {
        let file = File::create(path)?;
        let par_writer = ZBuilder::<Gzip, _>::new()
            .num_threads(threads)
            .from_writer(file);
            
        Ok(WriterType::MultiGz(Box::new(par_writer)))
    }
}