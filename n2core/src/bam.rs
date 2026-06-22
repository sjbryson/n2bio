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
            let byte: u8 = self.seq[i / 2];
            let nibble: u8 = if i % 2 == 0 {  // If the index is even, take the upper 4 bits. 
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
// BAM Flags trait and implementation
// ============================================================================

pub trait BamFlags {
    fn flag(&self) -> u16;
    
    // Bitwise default implementations
    #[inline] fn is_paired(&self)        -> bool { (self.flag() & 0x0001) != 0 }
    #[inline] fn is_proper(&self)        -> bool { (self.flag() & 0x0002) != 0 }
    #[inline] fn is_mapped(&self)        -> bool { (self.flag() & 0x0004) == 0 }
    #[inline] fn is_unmapped(&self)      -> bool { (self.flag() & 0x0004) != 0 }
    #[inline] fn is_mate_mapped(&self)   -> bool { (self.flag() & 0x0008) == 0 }
    #[inline] fn is_mate_unmapped(&self) -> bool { (self.flag() & 0x0008) != 0 }
    #[inline] fn is_revcomp(&self)       -> bool { (self.flag() & 0x0010) != 0 }
    #[inline] fn is_mate_revcomp(&self)  -> bool { (self.flag() & 0x0020) != 0 }
    #[inline] fn is_read1(&self)         -> bool { (self.flag() & 0x0040) != 0 }
    #[inline] fn is_read2(&self)         -> bool { (self.flag() & 0x0080) != 0 }
    #[inline] fn is_secondary(&self)     -> bool { (self.flag() & 0x0100) != 0 }
    #[inline] fn is_supplementary(&self) -> bool { (self.flag() & 0x0800) != 0 }
    // Derived logical methods
    #[inline] fn is_primary(&self) -> bool { 
        !self.is_secondary() && !self.is_supplementary() 
    }
}

impl BamFlags for BamRecord {
    #[inline]
    fn flag(&self) -> u16 {
        self.flag
    }
}

// ============================================================================
// BAM stats trait and implementation
// ============================================================================

pub trait BamStats {
    /// Retrieve an integer tag from the BAM tags byte array (e.g., b"AS" or b"NM")
    fn get_int_tag(&self, tag: &[u8; 2]) -> Option<i32>;
    
    /// AS (alignment score) / AL (alignment length)
    fn calculate_as_al(&self) -> Option<f32>;                 
    
    /// Sum of M, I, D, =, X operations from the CIGAR string
    fn calculate_alignment_length(&self) -> Option<u32>;      
    
    /// Alignment length / Read length
    fn calculate_alignment_proportion(&self) -> Option<f32>;  
    
    /// (Alignment length - NM) / Alignment length
    fn calculate_alignment_accuracy(&self) -> Option<f32>;    

    /// Calculates how many reference bases the alignment spans using the CIGAR string
    fn calculate_ref_span(&self) -> Option<u32>;
}

impl BamStats for BamRecord {
    
    fn get_int_tag(&self, tag: &[u8; 2]) -> Option<i32> {
        let mut i: usize = 0;
        let tags: &Vec<u8> = &self.tags;
        
        while i + 2 < tags.len() {
            let t1: u8 = tags[i];
            let t2: u8 = tags[i+1];
            let vtype: u8 = tags[i+2];
            i += 3;
            
            let is_target: bool = t1 == tag[0] && t2 == tag[1];
            
            match vtype {
                // 1-byte types: 'A' (char), 'c' (int8), 'C' (uint8)
                b'A' | b'c' | b'C' => {
                    if is_target {
                        return if vtype == b'c' {
                            Some(tags[i] as i8 as i32)
                        } else {
                            Some(tags[i] as i32)
                        };
                    }
                    i += 1;
                }
                // 2-byte types: 's' (int16), 'S' (uint16)
                b's' | b'S' => {
                    if i + 1 < tags.len() {
                        if is_target {
                            let bytes = [tags[i], tags[i+1]];
                            return if vtype == b's' {
                                Some(i16::from_le_bytes(bytes) as i32)
                            } else {
                                Some(u16::from_le_bytes(bytes) as i32)
                            };
                        }
                        i += 2;
                    } else { break; }
                }
                // 4-byte types: 'i' (int32), 'I' (uint32)
                b'i' | b'I' => {
                    if i + 3 < tags.len() {
                        if is_target {
                            let bytes: [u8; 4] = [tags[i], tags[i+1], tags[i+2], tags[i+3]];
                            return if vtype == b'i' {
                                Some(i32::from_le_bytes(bytes))
                            } else {
                                // Note: Converting u32 to i32. Safe for typical tags like AS/NM.
                                Some(u32::from_le_bytes(bytes) as i32) 
                            };
                        }
                        i += 4;
                    } else { break; }
                }
                // 4-byte float: 'f'
                b'f' => i += 4,
                // Null-terminated strings/hex: 'Z', 'H'
                b'Z' | b'H' => {
                    while i < tags.len() && tags[i] != 0 { i += 1; }
                    i += 1; // skip the null terminator
                }
                // Arrays: 'B' -> [type: 1 byte][count: 4 bytes][elements...]
                b'B' => {
                    if i + 4 < tags.len() {
                        let array_type: u8 = tags[i];
                        let count: usize = u32::from_le_bytes([tags[i+1], tags[i+2], tags[i+3], tags[i+4]]) as usize;
                        i += 5;
                        
                        let size: usize = match array_type {
                            b'c' | b'C' => 1,
                            b's' | b'S' => 2,
                            b'i' | b'I' | b'f' => 4,
                            _ => 0,
                        };
                        i += count * size;
                    } else { break; }
                }
                _ => break, // Unknown type, exit to prevent panics
            }
        }
        
        None
    }

    fn calculate_as_al(&self) -> Option<f32> {
        let align_len: u32 = self.calculate_alignment_length()?;
        let as_score: i32 = self.get_int_tag(b"AS")?;
        
        if align_len == 0 {
            None
        } else {
            Some(as_score as f32 / align_len as f32)
        }
    }

    fn calculate_alignment_length(&self) -> Option<u32> {
        let parsed_cigar: Vec<(CigarOp, u32)> = self.parsed_cigar();
        if parsed_cigar.is_empty() { return None; }
        
        let mut align_len: u32 = 0;
        for (op, len) in parsed_cigar {
            match op {
                // Operations that consume reference AND/OR mapping space
                CigarOp::Match | CigarOp::Insertion | CigarOp::Deletion | 
                CigarOp::SequenceMatch | CigarOp::SequenceMismatch => {
                    align_len += len;
                }
                _ => {} // Skip soft clips, hard clips, padding, etc.
            }
        }
        
        Some(align_len)
    }

    fn calculate_alignment_proportion(&self) -> Option<f32> {
        let align_len: u32 = self.calculate_alignment_length()?;
        
        // Use `l_seq` from the BAM fields for read length. 
        // If the sequence is missing ('*'), l_seq is 0, so use CIGAR.
        let read_len: u32 = if self.l_seq > 0 {
            self.l_seq as u32
        } else {
            let mut calculated_len: u32 = 0;
            for (op, len) in self.parsed_cigar() {
                match op {
                    CigarOp::Match | CigarOp::Insertion | CigarOp::SoftClip | 
                    CigarOp::SequenceMatch | CigarOp::SequenceMismatch => {
                        calculated_len += len;
                    }
                    _ => {}
                }
            }
            calculated_len
        };

        if read_len == 0 {
            None
        } else {
            Some(align_len as f32 / read_len as f32)
        }
    }

    fn calculate_alignment_accuracy(&self) -> Option<f32> {
        let align_len: u32 = self.calculate_alignment_length()?;
        let nm: i32 = self.get_int_tag(b"NM")?;
        
        if align_len == 0 {
            return None;
        }
        
        let matches: u32 = align_len.saturating_sub(nm as u32);
        Some((matches as f32 / align_len as f32) * 100.0)
    }

    fn calculate_ref_span(&self) -> Option<u32> {
        let parsed_cigar = self.parsed_cigar();
        if parsed_cigar.is_empty() { return None; }
        
        let span = parsed_cigar.iter()
            .filter(|(op, _)| matches!(
                op,
                CigarOp::Match | CigarOp::Deletion | CigarOp::Skip | 
                CigarOp::SequenceMatch | CigarOp::SequenceMismatch
            ))
            .map(|(_, len)| *len)
            .sum();
            
        Some(span)
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
    /// Smart open: Opens a BAM file utilizing the ReaderType buffer
    /// and handles the BGZF decompression block routing.
    pub fn open(path: &str) -> io::Result<Self> {
        let source: ReaderType = ReaderType::from_file(path)?;
        Ok(Self::new(source))
    }

    /// Initializes a new BamReader from a given ReaderType
    pub fn new(source: ReaderType) -> Self {
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

    /// Reads the next record into the `BamRecord` struct.
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
        let mut record: BamRecord = BamRecord::default();
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
        let mut record: BamRecord = BamRecord::default();
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

    /// Helper to create a BamRecord with specific CIGAR operations and tags
    fn create_record(l_seq: i32, cigar_ops: &[(u32, u8)], tags: Vec<u8>) -> BamRecord {
        let cigar: Vec<u32> = cigar_ops
            .iter()
            .map(|(len, op)| (len << 4) | (*op as u32))
            .collect();

        BamRecord {
            l_seq,
            cigar,
            tags,
            ..Default::default()
        }
    }

    #[test]
    fn test_get_int_tag_parsing() {
        // Tag format: [Tag1, Tag2, Type, Value...]
        // NM:C:5 (uint8), AS:i:1000 (int32), XX:Z:skip_me\0 (string)
        let tags: Vec<u8> = vec![
            b'X', b'X', b'Z', b's', b'k', b'i', b'p', 0, // String tag to skip
            b'N', b'M', b'C', 5,                         // uint8
            b'Z', b'Z', b'f', 0, 0, 0, 0,                // Float tag to skip
            b'A', b'S', b'i', 232, 3, 0, 0,              // int32 (1000 in little-endian)
            b'Y', b'Y', b'c', 254,                       // int8 (-2)
        ];

        let record: BamRecord = create_record(100, &[], tags);

        assert_eq!(record.get_int_tag(b"NM"), Some(5));
        assert_eq!(record.get_int_tag(b"AS"), Some(1000));
        assert_eq!(record.get_int_tag(b"YY"), Some(-2));
        assert_eq!(record.get_int_tag(b"XX"), None); // Not an int tag
        assert_eq!(record.get_int_tag(b"ZZ"), None); // Not an int tag
        assert_eq!(record.get_int_tag(b"NA"), None); // Missing tag
    }

    #[test]
    fn test_calculate_alignment_length() {
        // CIGAR: 50M (op 0), 10I (op 1), 5D (op 2), 20S (op 4)
        // Alignment length should be M + I + D = 50 + 10 + 5 = 65
        let record: BamRecord = create_record(0, &[(50, 0), (10, 1), (5, 2), (20, 4)], vec![]);
        
        assert_eq!(record.calculate_alignment_length(), Some(65));
    }

    #[test]
    fn test_calculate_as_al() {
        // CIGAR: 80M (op 0) -> align_len = 80
        // AS:s:120 (int16)
        let tags: Vec<u8> = vec![b'A', b'S', b's', 120, 0];
        let record: BamRecord = create_record(80, &[(80, 0)], tags);

        // 120.0 / 80.0 = 1.5
        assert_eq!(record.calculate_as_al(), Some(1.5));
    }

    #[test]
    fn test_calculate_alignment_proportion_with_l_seq() {
        // CIGAR: 50M (op 0), 50S (op 4) -> align_len = 50
        // l_seq = 100
        let record: BamRecord = create_record(100, &[(50, 0), (50, 4)], vec![]);

        // 50.0 / 100.0 = 0.5
        assert_eq!(record.calculate_alignment_proportion(), Some(0.5));
    }

    #[test]
    fn test_calculate_alignment_proportion_without_l_seq() {
        // Missing sequence data (l_seq = 0)
        // CIGAR: 50M (op 0), 10I (op 1), 20S (op 4) -> align_len = 60
        // Calculated read length: M + I + S = 50 + 10 + 20 = 80
        let record: BamRecord = create_record(0, &[(50, 0), (10, 1), (20, 4)], vec![]);

        // 60.0 / 80.0 = 0.75
        assert_eq!(record.calculate_alignment_proportion(), Some(0.75));
    }

    #[test]
    fn test_calculate_alignment_accuracy() {
        // CIGAR: 100M (op 0) -> align_len = 100
        // NM:C:5 (uint8) -> 5 mismatches
        let tags: Vec<u8> = vec![b'N', b'M', b'C', 5];
        let record: BamRecord = create_record(100, &[(100, 0)], tags);

        // (100 - 5) / 100 * 100 = 95.0%
        assert_eq!(record.calculate_alignment_accuracy(), Some(95.0));
    }

    #[test]
    fn test_calculate_alignment_accuracy_high_nm() {
        // NM is larger than alignment length (should safely saturate to 0 matches)
        // CIGAR: 10M -> align_len = 10
        // NM:C:15 -> 15 mismatches
        let tags: Vec<u8> = vec![b'N', b'M', b'C', 15];
        let record: BamRecord = create_record(10, &[(10, 0)], tags);

        assert_eq!(record.calculate_alignment_accuracy(), Some(0.0));
    }

    #[test]
    fn test_empty_cigar() {
        // Unmapped read with no CIGAR
        let record: BamRecord = create_record(100, &[], vec![b'A', b'S', b'C', 40]);

        assert_eq!(record.calculate_alignment_length(), None);
        assert_eq!(record.calculate_as_al(), None);
        assert_eq!(record.calculate_alignment_proportion(), None);
        assert_eq!(record.calculate_alignment_accuracy(), None);
    }

    #[test]
    fn test_calculate_ref_span() {
        // CIGAR: 50M (op 0), 10I (op 1), 5D (op 2), 100N (op 3), 20S (op 4)
        // Reference consuming operations are M (50), D (5), and N (100).
        // Total ref span should be 50 + 5 + 100 = 155
        // Insertions (I) and Soft Clips (S) do not consume the reference.
        let record = create_record(
            80, 
            &[(50, 0), (10, 1), (5, 2), (100, 3), (20, 4)], 
            vec![]
        );
        
        assert_eq!(record.calculate_ref_span(), Some(155));
    }
}