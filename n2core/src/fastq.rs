//! n2core/src/fastq.rs

use std::io::{self, Write, BufRead};
use std::marker::PhantomData;
use std::sync::Mutex;
use std::collections::HashMap;
use std::hash::{Hash, Hasher};
use rustc_hash::{FxHashMap, FxHasher};

use crate::sam::{SamRecord, SamStr, SamFlags};
use crate::readers::ReaderType;
use crate::writers::WriterType;
use crate::sequence::DnaSequence;

// ============================================================================
// FastqRecord <SingleRead> and <PairedRead> structs
// ============================================================================

pub trait FastqFormatter {
    fn format_id(id: &str) -> String;
}

#[derive(Debug, Clone)]
pub struct SingleRead;

impl FastqFormatter for SingleRead {
    fn format_id(id: &str) -> String { id.to_string() }
}

#[derive(Debug, Clone)]
pub struct Read1;

impl FastqFormatter for Read1 {
    fn format_id(id: &str) -> String { format!("{}/1", id) }
}

#[derive(Debug, Clone)]
pub struct Read2;

impl FastqFormatter for Read2 {
    fn format_id(id: &str) -> String { format!("{}/2", id) }
}

/// Struct representing an individual FASTQ read.
#[derive(Debug, Clone)]
pub struct FastqRecord<T: FastqFormatter> {
    pub id:   String,
    pub seq:  String,
    pub qual: String,
    _marker:  PhantomData<T>,
}

impl FastqRecord<SingleRead> {
    /// Create a single FASTQ record from a SamStr.
    pub fn from_samstr(sam: &SamStr) -> Self {
        let (qname, seq, qual) = sam.to_fastq_fields();
        Self::new(qname.to_string(), seq.to_string(), qual.to_string())
    }
    /// Create a single FASTQ record from a SamRecord.
    pub fn from_samrec(rec: &SamRecord) -> Self {
        Self::new(rec.qname.clone(), rec.seq.clone(), rec.qual.clone())
    }
}

impl<T: FastqFormatter> FastqRecord<T> {
    pub fn new(id: String, seq: String, qual: String) -> Self {
        Self { id, seq, qual, _marker: PhantomData }
    }
}

/// Struct for a paired FASTQ read
#[derive(Debug, Clone)]
pub struct PairedFastqRecord {
    pub r1: FastqRecord<Read1>,
    pub r2: FastqRecord<Read2>,
}

/// Enum to parse reads to correct orientation.
pub enum PairedRead {
    R1(FastqRecord<Read1>),
    R2(FastqRecord<Read2>),
}

impl PairedRead {
    pub fn id(&self) -> &str {
        match self {
            PairedRead::R1(rec) => &rec.id,
            PairedRead::R2(rec) => &rec.id,
        }
    }

    pub fn from_samstr(sam: &SamStr) -> Self {
        let (qname, seq, qual) = sam.to_fastq_fields();

        if sam.is_read1() {
            PairedRead::R1(FastqRecord::<Read1>::new(
                qname.to_string(), seq.to_string(), qual.to_string()
            ))
        } else {
            PairedRead::R2(FastqRecord::<Read2>::new(
                qname.to_string(), seq.to_string(), qual.to_string()
            ))
        }
    }

    pub fn from_samrec(mut rec: SamRecord) -> Self {
        // Reverse-complement if the 0x10 flag is set
        if rec.is_revcomp() {
            rec.seq  = rec.seq.reverse_complement();
            rec.qual = rec.qual.chars().rev().collect();
        }

        if rec.is_read1() {
            PairedRead::R1(FastqRecord::<Read1>::new(
                rec.qname, rec.seq, rec.qual 
            ))
        } else {
            PairedRead::R2(FastqRecord::<Read2>::new(
                rec.qname, rec.seq, rec.qual
            ))
        }
    }
}


// ============================================================================
// Writers
// ============================================================================

/// A writer for standard single-end FASTQ records.
pub struct FastqWriter<W: Write> {
    pub writer: W,
}

impl<W: Write> FastqWriter<W> {
    pub fn new(writer: W) -> Self {
        Self { writer }
    }

    /// Writes a record using raw bytes to bypass `std::fmt` overhead.
    pub fn write_record<T: FastqFormatter>(&mut self, rec: &FastqRecord<T>) -> io::Result<()> {
        self.writer.write_all(b"@")?;
        self.writer.write_all(T::format_id(&rec.id).as_bytes())?;
        self.writer.write_all(b"\n")?;
        self.writer.write_all(rec.seq.as_bytes())?;
        self.writer.write_all(b"\n+\n")?;
        self.writer.write_all(rec.qual.as_bytes())?;
        self.writer.write_all(b"\n")
    }

    /// Flushes the underlying writer.
    pub fn flush(&mut self) -> io::Result<()> {
        self.writer.flush()
    }

    /// Consumes the FastqWriter, returning the underlying writer.
    pub fn into_inner(self) -> W {
        self.writer
    }
}

impl FastqWriter<WriterType> {
    /// Automatically detects file type based on extension.
    pub fn create(path: &str) -> io::Result<Self> {
        Ok(Self::new(WriterType::create(path)?))
    }

    /// Uses multi-threaded gzip for `.gz` files.
    pub fn create_with_threads(path: &str, threads: usize) -> io::Result<Self> {
        let writer: WriterType = if path.ends_with(".gz") {
            WriterType::to_multithreaded_gz(path, threads)?
        } else {
            WriterType::create(path)?
        };
        Ok(Self::new(writer))
    }
}

/// A writer for paired-end FASTQ records, managing two output streams simultaneously.
pub struct PairedFastqWriter<W1: Write, W2: Write> {
    pub w1: FastqWriter<W1>,
    pub w2: FastqWriter<W2>,
}

impl<W1: Write, W2: Write> PairedFastqWriter<W1, W2> {
    pub fn new(writer1: W1, writer2: W2) -> Self {
        Self {
            w1: FastqWriter::new(writer1),
            w2: FastqWriter::new(writer2),
        }
    }

    pub fn write_pair(&mut self, pair: &PairedFastqRecord) -> io::Result<()> {
        self.w1.write_record(&pair.r1)?;
        self.w2.write_record(&pair.r2)?;
        Ok(())
    }

    /// Flushes both R1 and R2 underlying writers.
    pub fn flush(&mut self) -> io::Result<()> {
        self.w1.flush()?;
        self.w2.flush()
    }

    /// Consumes the PairedFastqWriter, returning the two underlying writers.
    pub fn into_inner(self) -> (W1, W2) {
        (self.w1.into_inner(), self.w2.into_inner())
    }
}

impl PairedFastqWriter<WriterType, WriterType> {
    /// Automatically detects file type based on extension.
    pub fn create(path1: &str, path2: &str) -> io::Result<Self> {
        Ok(Self::new(WriterType::create(path1)?, WriterType::create(path2)?))
    }

    /// Uses multi-threaded gzip for `.gz` files.
    pub fn create_with_threads(path1: &str, path2: &str, threads: usize) -> io::Result<Self> {
        let w1: WriterType = if path1.ends_with(".gz") { WriterType::to_multithreaded_gz(path1, threads)? } else { WriterType::create(path1)? };
        let w2: WriterType = if path2.ends_with(".gz") { WriterType::to_multithreaded_gz(path2, threads)? } else { WriterType::create(path2)? };
        Ok(Self::new(w1, w2))
    }
}

/// A writer for interleaved paired-end FASTQ records.
pub struct InterleavedFastqWriter<W: Write> {
    pub writer: FastqWriter<W>,
}

impl<W: Write> InterleavedFastqWriter<W> {
    pub fn new(writer: W) -> Self {
        Self {
            writer: FastqWriter::new(writer),
        }
    }

    /// Writes R1 then R2 to the same stream.
    pub fn write_pair(&mut self, pair: &PairedFastqRecord) -> io::Result<()> {
        self.writer.write_record(&pair.r1)?;
        self.writer.write_record(&pair.r2)?;
        Ok(())
    }

    /// Flushes the underlying writer.
    pub fn flush(&mut self) -> io::Result<()> {
        self.writer.flush()
    }

    /// Consumes the InterleavedFastqWriter, returning the underlying writer.
    pub fn into_inner(self) -> W {
        self.writer.into_inner()
    }
}

impl InterleavedFastqWriter<WriterType> {
    pub fn create(path: &str) -> io::Result<Self> {
        Ok(Self::new(WriterType::create(path)?))
    }

    pub fn create_with_threads(path: &str, threads: usize) -> io::Result<Self> {
        let writer: WriterType = if path.ends_with(".gz") { WriterType::to_multithreaded_gz(path, threads)? } else { WriterType::create(path)? };
        Ok(Self::new(writer))
    }
}

// ============================================================================
// Readers
// ============================================================================

/// A reader for parsing individual FASTQ records.
pub struct FastqReader<R: BufRead, T> {
    pub reader: R,
    _marker: PhantomData<T>,
}

impl<R: BufRead, T> FastqReader<R, T> {
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            _marker: PhantomData,
        }
    }
}

impl<T> FastqReader<ReaderType, T>{
    /// Automatically detects file type based on extension.
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

    pub fn from_bz(path: &str) -> io::Result<Self> {
        Ok(Self::new(ReaderType::from_bz(path)?))
    }
}

impl<R: BufRead, T: FastqFormatter> Iterator for FastqReader<R, T> {
    type Item = io::Result<FastqRecord<T>>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut id_line = String::new();
        
        // 1. Read the ID line.
        match self.reader.read_line(&mut id_line) {
            Ok(0) => return None, // Clean EOF
            Ok(_) => {},
            Err(e) => return Some(Err(e)),
        }

        let mut seq:  String = String::new();
        let mut plus: String = String::new();
        let mut qual: String = String::new();

        // 2. Read the next 3 lines. If any fail, it's an unexpected EOF/truncation.
        if let Err(e) = self.reader.read_line(&mut seq) { return Some(Err(e)); }
        if let Err(e) = self.reader.read_line(&mut plus) { return Some(Err(e)); }
        if let Err(e) = self.reader.read_line(&mut qual) { return Some(Err(e)); }

        // 3. Clean up the strings.
        let id:   String = id_line.trim_start_matches('@').trim_end().to_string();
        let seq:  String = seq.trim_end().to_string();
        let qual: String = qual.trim_end().to_string();

        // 4. Validate that the record wasn't truncated in the middle of a file.
        if seq.is_empty() || plus.is_empty() || qual.is_empty() {
            return Some(Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "Truncated FASTQ record: missing sequence or quality lines",
            )));
        }

        // 5. Return the typed FastqRecord.
        Some(Ok(FastqRecord::new(id, seq, qual)))
    }
}

/// A reader that steps through two FASTQ files simultaneously, yielding paired records.
pub struct PairedFastqReader<R1: BufRead, R2: BufRead> {
    pub r1_reader: FastqReader<R1, Read1>,
    pub r2_reader: FastqReader<R2, Read2>,
}

impl<R1: BufRead, R2: BufRead> PairedFastqReader<R1, R2> {
    pub fn new(reader1: R1, reader2: R2) -> Self {
        Self {
            r1_reader: FastqReader::new(reader1),
            r2_reader: FastqReader::new(reader2),
        }
    }
}

impl PairedFastqReader<ReaderType, ReaderType> {
    /// Opens two auto-detecting FASTQ streams.
    pub fn open(path1: &str, path2: &str) -> io::Result<Self> {
        Ok(Self::new(ReaderType::open(path1)?, ReaderType::open(path2)?))
    }

    /// Opens two uncompressed FASTQ files and returns a paired reader.
    pub fn from_files(path1: &str, path2: &str) -> io::Result<Self> {
        let r1: ReaderType = ReaderType::from_file(path1)?;
        let r2: ReaderType = ReaderType::from_file(path2)?;
        Ok(Self::new(r1, r2))
    }

    /// Opens two gzipped FASTQ files and returns a paired reader.
    pub fn from_gzs(path1: &str, path2: &str) -> io::Result<Self> {
        let r1: ReaderType = ReaderType::from_gz(path1)?;
        let r2: ReaderType = ReaderType::from_gz(path2)?;
        Ok(Self::new(r1, r2))
    }

    /// Opens two bz2 compressed FASTQ files and returns a paired reader.
    pub fn from_bzs(path1: &str, path2: &str) -> io::Result<Self> {
        let r1: ReaderType = ReaderType::from_bz(path1)?;
        let r2: ReaderType = ReaderType::from_bz(path2)?;
        Ok(Self::new(r1, r2))
    }
}

impl<R1: BufRead, R2: BufRead> Iterator for PairedFastqReader<R1, R2> {
    type Item = io::Result<PairedFastqRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        let rec1: Option<Result<FastqRecord<Read1>, io::Error>> = self.r1_reader.next();
        let rec2: Option<Result<FastqRecord<Read2>, io::Error>> = self.r2_reader.next();

        match (rec1, rec2) {
            // Both read successfully
            (Some(Ok(r1)), Some(Ok(r2))) => Some(Ok(PairedFastqRecord { r1, r2 })),
            
            // Both hit clean EOF simultaneously
            (None, None) => None, 
            
            // Standard I/O Errors bubble up
            (Some(Err(e)), _) | (_, Some(Err(e))) => Some(Err(e)), 
            
            // One file ended before the other - orphaned pairs
            _ => Some(Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Mismatched FASTQ files: one file ended before the other",
            ))),
        }
    }
}

/// A reader that parses interleaved paired-end FASTQ records from a single stream.
pub struct InterleavedFastqReader<R: BufRead> {
    pub reader: FastqReader<R, SingleRead>,
}

impl<R: BufRead> InterleavedFastqReader<R> {
    pub fn new(reader: R) -> Self {
        Self {
            reader: FastqReader::new(reader),
        }
    }
}

// Convenience constructors using ReaderType
impl InterleavedFastqReader<ReaderType> {
    /// Automatically detects file type based on extension.
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

    pub fn from_bz(path: &str) -> io::Result<Self> {
        Ok(Self::new(ReaderType::from_bz(path)?))
    }
}

impl<R: BufRead> Iterator for InterleavedFastqReader<R> {
    type Item = io::Result<PairedFastqRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        // Read the first record (R1)
        let rec1: FastqRecord<SingleRead> = match self.reader.next() {
            Some(Ok(r)) => r,
            Some(Err(e)) => return Some(Err(e)),
            None => return None, // Clean EOF at the start of a pair
        };

        // Read the second record (R2)
        let rec2: FastqRecord<SingleRead> = match self.reader.next() {
            Some(Ok(r)) => r,
            Some(Err(e)) => return Some(Err(e)),
            None => {
                return Some(Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Truncated interleaved FASTQ: Found R1 without a matching R2",
                )))
            }
        };

        // Convert SingleReads to Read1 and Read2 to create PairedFastqRecord struct
        let r1: FastqRecord<Read1> = FastqRecord::<Read1>::new(rec1.id, rec1.seq, rec1.qual);
        let r2: FastqRecord<Read2> = FastqRecord::<Read2>::new(rec2.id, rec2.seq, rec2.qual);

        Some(Ok(PairedFastqRecord { r1, r2 }))
    }
}

// ============================================================================
// Mate pairing HashMaps
// ============================================================================

/// MateMap - a struct for using a HashMap to collect mate pairs
pub struct MateMap {
    pending: HashMap<String, PairedRead>,
}

impl MateMap {
    pub fn new() -> Self { Self { pending: HashMap::new() } }

    pub fn process(&mut self, rec: PairedRead) -> Option<PairedFastqRecord> {
        let qname: String = rec.id().to_string();

        if let Some(mate) = self.pending.remove(&qname) {
            let pair: PairedFastqRecord = match (rec, mate) {
                (PairedRead::R1(r1), PairedRead::R2(r2)) |
                (PairedRead::R2(r2), PairedRead::R1(r1)) => {
                    PairedFastqRecord { r1, r2 }
                },
                _ => {
                    eprintln!("Warning: Malformed SAM pairs for QNAME {}", qname);
                    return None; 
                }
            };
            Some(pair)
        } else {
            self.pending.insert(qname, rec);
            None
        }
    }

    pub fn orphan_count(&self) -> usize {
        self.pending.len()
    }
}

/// A highly concurrent, sharded map for collecting mate pairs across multiple threads.
pub struct ShardedMateMap {
    shards: Vec<Mutex<FxHashMap<String, PairedRead>>>,
    num_shards: usize,
}

impl ShardedMateMap {
    /// Creates a new ShardedMateMap. 
    /// A good rule of thumb is to set num_shards to at least 4x to 8x 
    /// the number of threads you plan to run. If you are using a 16-core machine, 
    /// 64 or 128 shards virtually guarantee that two threads will almost never 
    /// try to lock the exact same shard at the same time.
    pub fn new(num_shards: usize) -> Self {
        let mut shards: Vec<Mutex<HashMap<String, PairedRead, std::hash::BuildHasherDefault<FxHasher>>>> = Vec::with_capacity(num_shards);
        for _ in 0..num_shards {
            // FxHashMap::default() uses the fast FxHasher automatically
            shards.push(Mutex::new(FxHashMap::default()));
        }
        
        Self { shards, num_shards }
    }

    /// Determines which shard a given QNAME belongs to.
    #[inline]
    fn get_shard_index(&self, qname: &str) -> usize {
        let mut hasher: FxHasher = FxHasher::default();
        qname.hash(&mut hasher);
        (hasher.finish() as usize) % self.num_shards
    }

    /// Processes a read safely across threads.
    pub fn process(&self, rec: PairedRead) -> Option<PairedFastqRecord> {
        let shard_idx: usize = self.get_shard_index(rec.id());

        // 1. Lock the specific shard.
        let mut shard: std::sync::MutexGuard<'_, HashMap<String, PairedRead, std::hash::BuildHasherDefault<FxHasher>>> = self.shards[shard_idx].lock().unwrap();

        // 2. Look up using the reference.
        if let Some(mate) = shard.remove(rec.id()) {
            let pair: PairedFastqRecord = match (rec, mate) {
                (PairedRead::R1(r1), PairedRead::R2(r2)) |
                (PairedRead::R2(r2), PairedRead::R1(r1)) => {
                    PairedFastqRecord { r1, r2 }
                },
                (bad_rec, _bad_mate) => {
                    // Bind the moved record to bad_rec to print its ID
                    eprintln!("Warning: Malformed SAM pairs for QNAME {}", bad_rec.id());
                    return None;
                }
            };
            Some(pair)
        } else {
            // 3. Only allocate a String to insert it.
            let qname_owned: String = rec.id().to_string();
            shard.insert(qname_owned, rec);
            None
        }
    }

    /// Safely aggregates the total orphan count across all shards.
    pub fn orphan_count(&self) -> usize {
        self.shards
            .iter()
            .map(|shard| shard.lock().unwrap().len())
            .sum()
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_formatters() {
        assert_eq!(SingleRead::format_id("read1"), "read1");
        assert_eq!(Read1::format_id("read1"), "read1/1");
        assert_eq!(Read2::format_id("read1"), "read1/2");
    }

    #[test]
    fn test_fastq_reader() {
        // Mock a simple FASTQ file in memory
        let fq_data = b"@read_A\nACGT\n+\n!!!!\n@read_B\nTGCA\n+\n####\n";
        let cursor = Cursor::new(fq_data);
        let mut reader = FastqReader::<_, SingleRead>::new(cursor);
        
        let rec1 = reader.next().unwrap().unwrap();
        assert_eq!(rec1.id, "read_A");
        assert_eq!(rec1.seq, "ACGT");
        assert_eq!(rec1.qual, "!!!!");

        let rec2 = reader.next().unwrap().unwrap();
        assert_eq!(rec2.id, "read_B");
        assert_eq!(rec2.seq, "TGCA");
        assert_eq!(rec2.qual, "####");
        
        assert!(reader.next().is_none());
    }

    #[test]
    fn test_truncated_fastq() {
        let fq_data = b"@read_A\nACGT\n+\n"; // Missing qual line
        let cursor = Cursor::new(fq_data);
        let mut reader = FastqReader::<_, SingleRead>::new(cursor);
        
        let result = reader.next().unwrap();
        assert!(result.is_err());
        assert_eq!(result.unwrap_err().kind(), io::ErrorKind::UnexpectedEof);
    }

    #[test]
    fn test_mate_map() {
        let mut map = MateMap::new();
        let r1 = PairedRead::R1(FastqRecord::<Read1>::new("pair1".to_string(), "A".to_string(), "!".to_string()));
        let r2 = PairedRead::R2(FastqRecord::<Read2>::new("pair1".to_string(), "C".to_string(), "!".to_string()));
        
        // Processing R1 first should cache it and return None
        assert!(map.process(r1).is_none());
        assert_eq!(map.orphan_count(), 1);
        
        // Processing R2 should find the match and return the pair
        let pair = map.process(r2).unwrap();
        assert_eq!(pair.r1.seq, "A");
        assert_eq!(pair.r2.seq, "C");
        assert_eq!(map.orphan_count(), 0);
    }
}