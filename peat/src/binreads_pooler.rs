//! n2bio/peat/src/binpooler.rs
//! 

use std::collections::HashMap;
use std::io;
use std::num::NonZeroUsize;
use std::path::PathBuf;
use lru::LruCache;

// Importing your structs
use n2bio::fastq::{PairedFastqRecord, PairedFastqWriter};
use n2bio::writers::WriterType;

pub struct BinnedFastqPool {
    /// LRU Cache holding open active PairedFastqWriters per bin.
    cache: LruCache<String, PairedFastqWriter<WriterType, WriterType>>,
    
    /// In-memory record buffers per bin to batch disk writes.
    buffers: HashMap<String, Vec<PairedFastqRecord>>,
    
    /// Base output directory for generated bin files.
    out_dir: PathBuf,
    
    /// Maximum pairs to hold in memory per bin before flushing to disk.
    buffer_capacity: usize,
    
    /// Number of threads dedicated to gzip compression engine per writer stream.
    gz_threads: usize,
    
    /// Custom file prefix or naming template.
    prefix: String,
}

impl BinnedFastqPool {
    pub fn new(
        max_open_bins: usize,
        buffer_capacity: usize,
        out_dir: PathBuf,
        gz_threads: usize,
        prefix: String,
    ) -> Self {
        Self {
            cache: LruCache::new(NonZeroUsize::new(max_open_bins).expect("max_open_bins must be > 0")),
            buffers: HashMap::new(),
            out_dir,
            buffer_capacity,
            gz_threads,
            prefix,
        }
    }

    /// Pushes a completed paired read record into the specified bin's buffer.
    /// Flushes to disk if the bin's buffer capacity threshold is reached.
    pub fn push(&mut self, bin_id: &str, pair: PairedFastqRecord) -> io::Result<()> {
        let buf: &mut Vec<PairedFastqRecord> = self.buffers.entry(bin_id.to_string()).or_insert_with(Vec::new);
        buf.push(pair);

        if buf.len() >= self.buffer_capacity {
            self.flush_bin(bin_id)?;
        }
        Ok(())
    }

    /// Helper to flush a specific bin's memory buffer into its corresponding `PairedFastqWriter`.
    fn flush_bin(&mut self, bin_id: &str) -> io::Result<()> {
        if let Some(reads) = self.buffers.get_mut(bin_id) {
            if reads.is_empty() {
                return Ok(());
            }

            // If the writer isn't in the active LRU cache, instantiate it
            if !self.cache.contains(bin_id) {
                let r1_path: PathBuf = self.out_dir.join(format!("{}.{}.r1.fq.gz", self.prefix, bin_id));
                let r2_path: PathBuf = self.out_dir.join(format!("{}.{}.r2.fq.gz", self.prefix, bin_id));

                let writer: PairedFastqWriter<WriterType, WriterType> = PairedFastqWriter::create_with_threads(
                    r1_path.to_str().expect("Invalid UTF-8 in output path"),
                    r2_path.to_str().expect("Invalid UTF-8 in output path"),
                    self.gz_threads,
                )?;

                // Push to cache. If an old writer is evicted flush its buffers.
                if let Some((_evicted_id, mut evicted_writer)) = self.cache.push(bin_id.to_string(), writer) {
                    evicted_writer.flush()?;
                }
            }

            // Get writer from cache and drain memory buffer to disk
            let writer: &mut PairedFastqWriter<WriterType, WriterType> = self.cache.get_mut(bin_id).unwrap();
            for pair in reads.drain(..) {
                writer.write_pair(&pair)?;
            }
        }
        Ok(())
    }

    /// Finalize the pool: flush remaining memory buffers and all active stream writers in the LRU cache.
    pub fn finish_all(&mut self) -> io::Result<()> {
        // Flush buffered reads in memory
        let bin_ids: Vec<String> = self.buffers.keys().cloned().collect();
        for bin_id in bin_ids {
            self.flush_bin(&bin_id)?;
        }

        // Flush remaining writers in the LRU cache
        for (_bin_id, writer) in self.cache.iter_mut() {
            writer.flush()?;
        }

        Ok(())
    }
}