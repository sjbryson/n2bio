//! n2bio/n2core/src/bam.rs

use byteorder::{LittleEndian, ReadBytesExt};
use std::io::{self, Read, ErrorKind};
use flate2::read::MultiGzDecoder;
use crate::readers::ReaderType;

// BAM 16-character sequence alphabet (Index matches the 4-bit integer).
// 0:=, 1:A, 2:C, 3:M, 4:G, 5:R, 6:S, 7:V, 8:T, 9:W, 10:Y, 11:H, 12:K, 13:D, 14:B, 15:N
const BAM_SEQ_LOOKUP: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

// ============================================================================
// Structs & Enums
// ============================================================================

#[derive(Debug, Clone)]
pub struct BamHeader {
    pub text: String,                       // Plain-text SAM header.
    pub references: Vec<ReferenceSequence>, // Reference sequence dictionary
}

#[derive(Debug, Clone)]
pub struct ReferenceSequence {
    pub name: String,  // Reference name.
    pub length: i32,   // Reference length.
}

#[derive(Debug, PartialEq, Eq)]
pub enum CigarOp {
    Match,             // 0: M
    Insertion,         // 1: I
    Deletion,          // 2: D
    Skip,              // 3: N
    SoftClip,          // 4: S
    HardClip,          // 5: H
    Padding,           // 6: P
    SequenceMatch,     // 7: =
    SequenceMismatch,  // 8: X
    Unknown(u8),       // Catch-all for unexpected values
}

// ============================================================================
// BAM Record struct and implementation
// ============================================================================

/// Represents a single BAM alignment record.
/// The fields correspond directly to the 32-byte fixed-length block 
/// and the subsequent variable-length data in the BAM spec.
#[derive(Debug, Default, Clone)]
pub struct BamRecord {
    // Fixed Length Fields  //
    pub ref_id: i32,        // Reference sequence ID (-1 for unmapped)
    pub pos: i32,           // 0-based leftmost coordinate
    pub l_read_name: u8,    // Length of the read name (including null byte)
    pub mapq: u8,           // Mapping quality
    pub bin: u16,           // BAM index bin
    pub n_cigar_op: u16,    // Num CIGAR operations
    pub flag: u16,          // Bitwise flags
    pub l_seq: i32,         // Seq length
    pub next_ref_id: i32,   // Reference ID for next segment (mate)
    pub next_pos: i32,      // 0-based leftmost coordinate of the next segment
    pub tlen: i32,          // Template length (insert size)
    // Variable Length Data //
    pub read_name: Vec<u8>, // Read name (null-terminated ASCII)
    pub cigar: Vec<u32>,    // CIGAR operations encoded as 32-bit int
    pub seq: Vec<u8>,       // 4-bit bases (e.g., =0, A=1, C=2, G=4, T=8, N=15)
    pub qual: Vec<u8>,      // Phred quality (qual + 33, but 0-indexed here)
    pub tags: Vec<u8>,      // Tags like NM:i:1, MD:Z:...)
}

impl BamRecord {
    /// Decodes the 4-bit packed sequence into standard ASCII bytes (e.g., b"ACGT")
    pub fn decoded_sequence(&self) -> Vec<u8> {
        let mut ascii_seq: Vec<u8> = Vec::with_capacity(self.l_seq as usize);

        for i in 0..(self.l_seq as usize) {
            let byte = self.seq[i / 2];
            let nibble = if i % 2 == 0 {  // If the index is even, take the upper 4 bits. 
                (byte >> 4) & 0x0F
            } else {
                byte & 0x0F                   // If odd, take the lower 4 bits.
            };

            ascii_seq.push(BAM_SEQ_LOOKUP[nibble as usize]);
        }

        ascii_seq
    }
    
    /// Unpacks binary CIGAR array into a Vector of (Operation, Length) tuples
    pub fn parsed_cigar(&self) -> Vec<(CigarOp, u32)> {
        let mut parsed: Vec<(CigarOp, u32)> = Vec::with_capacity(self.n_cigar_op as usize);

        for &op_data in &self.cigar {
            let length: u32 = op_data >> 4;           // Shift the 4-bit operation.
            let op_code: u8 = (op_data & 0x0F) as u8; // Mask all but last 4 bits.

            let op = match op_code {
                0 => CigarOp::Match,
                1 => CigarOp::Insertion,
                2 => CigarOp::Deletion,
                3 => CigarOp::Skip,
                4 => CigarOp::SoftClip,
                5 => CigarOp::HardClip,
                6 => CigarOp::Padding,
                7 => CigarOp::SequenceMatch,
                8 => CigarOp::SequenceMismatch,
                _ => CigarOp::Unknown(op_code),
            };

            parsed.push((op, length));
        }

        parsed
    }
}

// ============================================================================
// BAM Reader
// ============================================================================

/// BAM Reader that handles BGZF decompression
pub struct BamReader {
    // Wrap input in BufReader and pass to MultiGzDecoder to handle BGZF blocks.
    bgzf_stream: MultiGzDecoder<ReaderType>,
}

impl BamReader {
    /// Smart open: Opens a BAM file utilizing the optimized ReaderType buffer
    /// and seamlessly handles the BGZF decompression block routing.
    pub fn open(path: &str) -> io::Result<Self> {
        let source: ReaderType = ReaderType::from_file(path)?;
        Ok(Self::new(source))
    }

    /// Initializes a new BamReader from a given ReaderType
    pub fn new(source: ReaderType) -> Self {
        // ReaderType is already heavily buffered; we pass it directly into the decoder.
        let bgzf_stream: MultiGzDecoder<ReaderType> = MultiGzDecoder::new(source);
        Self { bgzf_stream }
    }

    /// Reads and parses the BAM header and reference dictionary
    pub fn read_header(&mut self) -> io::Result<BamHeader> {
        let mut magic: [u8; 4] = [0u8; 4];
        self.bgzf_stream.read_exact(&mut magic)?;
        if &magic != b"BAM\x01" {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Invalid BAM magic string: File is not a valid BAM",
            ));
        }

        let l_text: i32 = self.bgzf_stream.read_i32::<LittleEndian>()?;
        let mut text_bytes: Vec<u8> = vec![0u8; l_text as usize];
        self.bgzf_stream.read_exact(&mut text_bytes)?;
        
        let text: String = String::from_utf8(text_bytes).map_err(|e| {
            io::Error::new(io::ErrorKind::InvalidData, format!("Invalid UTF-8 in header: {}", e))
        })?;

        let n_ref: i32 = self.bgzf_stream.read_i32::<LittleEndian>()?;
        let mut references: Vec<ReferenceSequence> = Vec::with_capacity(n_ref as usize);

        for _ in 0..n_ref {
            let l_name: i32 = self.bgzf_stream.read_i32::<LittleEndian>()?;
            let mut name_bytes: Vec<u8> = vec![0u8; l_name as usize];
            self.bgzf_stream.read_exact(&mut name_bytes)?;
            
            if let Some(b'\0') = name_bytes.last() {
                name_bytes.pop();
            }
            
            let name: String = String::from_utf8(name_bytes).map_err(|e| {
                io::Error::new(io::ErrorKind::InvalidData, format!("Invalid UTF-8 in ref name: {}", e))
            })?;

            let l_ref: i32 = self.bgzf_stream.read_i32::<LittleEndian>()?;

            references.push(ReferenceSequence {
                name,
                length: l_ref,
            });
        }

        Ok(BamHeader { text, references })
    }

    /// Reads the next record into the provided `BamRecord` struct.
    /// Returns `Ok(true)` if a record was read, and `Ok(false)` on EOF.
    pub fn read_record(&mut self, record: &mut BamRecord) -> io::Result<bool> {
        let block_size: i32 = match self.bgzf_stream.read_i32::<LittleEndian>() {
            Ok(size) => size,
            Err(ref e) if e.kind() == ErrorKind::UnexpectedEof => return Ok(false), 
            Err(e) => return Err(e),
        };

        record.ref_id      = self.bgzf_stream.read_i32::<LittleEndian>()?;
        record.pos         = self.bgzf_stream.read_i32::<LittleEndian>()?;
        record.l_read_name = self.bgzf_stream.read_u8()?;
        record.mapq        = self.bgzf_stream.read_u8()?;
        record.bin         = self.bgzf_stream.read_u16::<LittleEndian>()?;
        record.n_cigar_op  = self.bgzf_stream.read_u16::<LittleEndian>()?;
        record.flag        = self.bgzf_stream.read_u16::<LittleEndian>()?;
        record.l_seq       = self.bgzf_stream.read_i32::<LittleEndian>()?;
        record.next_ref_id = self.bgzf_stream.read_i32::<LittleEndian>()?;
        record.next_pos    = self.bgzf_stream.read_i32::<LittleEndian>()?;
        record.tlen        = self.bgzf_stream.read_i32::<LittleEndian>()?;

        record.read_name.clear();
        record.cigar.clear();
        record.seq.clear();
        record.qual.clear();
        record.tags.clear();

        record.read_name.resize(record.l_read_name as usize, 0);
        self.bgzf_stream.read_exact(&mut record.read_name)?;

        record.cigar.reserve(record.n_cigar_op as usize);
        for _ in 0..record.n_cigar_op {
            record.cigar.push(self.bgzf_stream.read_u32::<LittleEndian>()?);
        }

        let seq_bytes_len: usize = ((record.l_seq + 1) / 2) as usize;
        record.seq.resize(seq_bytes_len, 0);
        self.bgzf_stream.read_exact(&mut record.seq)?;

        record.qual.resize(record.l_seq as usize, 0);
        self.bgzf_stream.read_exact(&mut record.qual)?;

        let cigar_bytes_len: i32 = (record.n_cigar_op as i32) * 4;
        let bytes_read_so_far: i32 = 32 + (record.l_read_name as i32) + cigar_bytes_len + (seq_bytes_len as i32) + record.l_seq;
        let tags_len: i32 = block_size - bytes_read_so_far;
        
        if tags_len > 0 {
            record.tags.resize(tags_len as usize, 0);
            self.bgzf_stream.read_exact(&mut record.tags)?;
        }

        Ok(true)
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decoded_sequence_even_length() {
        let mut record = BamRecord::default();
        // Construct the expected sequence "ACGT"
        // A=1, C=2  => (1 << 4) | 2 = 18 (0x12)
        // G=4, T=8  => (4 << 4) | 8 = 72 (0x48)
        record.l_seq = 4;
        record.seq = vec![0x12, 0x48];
        
        let decoded: Vec<u8> = record.decoded_sequence();
        assert_eq!(decoded, b"ACGT");
    }

    #[test]
    fn test_decoded_sequence_odd_length() {
        let mut record = BamRecord::default();
        // Construct the expected sequence "ACG"
        // A=1, C=2  => (1 << 4) | 2 = 18 (0x12)
        // G=4, N=15 => (4 << 4) | 15 = 79 (0x4F) -> The 'N' is a padded ignored byte
        record.l_seq = 3;
        record.seq = vec![0x12, 0x4F];
        
        let decoded: Vec<u8> = record.decoded_sequence();
        assert_eq!(decoded, b"ACG");
    }

    #[test]
    fn test_parsed_cigar() {
        let mut record: BamRecord = BamRecord::default();
        record.n_cigar_op = 3;
        
        // 1. 3M (3 Matches)     => Op 0, Len 3 -> (3 << 4) | 0 = 48  (0x30)
        // 2. 1I (1 Insertion)   => Op 1, Len 1 -> (1 << 4) | 1 = 17  (0x11)
        // 3. 2D (2 Deletions)   => Op 2, Len 2 -> (2 << 4) | 2 = 34  (0x22)
        record.cigar = vec![0x30, 0x11, 0x22];

        let parsed: Vec<(CigarOp, u32)> = record.parsed_cigar();
        
        assert_eq!(parsed.len(), 3);
        assert_eq!(parsed[0], (CigarOp::Match, 3));
        assert_eq!(parsed[1], (CigarOp::Insertion, 1));
        assert_eq!(parsed[2], (CigarOp::Deletion, 2));
    }
}