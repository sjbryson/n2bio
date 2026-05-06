//! n2core/src/sam.rs
//! Core structures and traits for parsing, manipulating, and analyzing SAM records.

use std::borrow::Cow;
use std::convert::TryFrom;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

use crate::fastq::{FastqFormatter, FastqRecord};
use crate::readers::ReaderType;
use crate::sequence::DnaSequence;
use thiserror::Error;

// ============================================================================
// Errors
// ============================================================================

/// Custom errors for SAM processing
#[derive(Error, Debug)]
pub enum SamError {
    #[error("IO Error: {0}")]
    Io(#[from] io::Error),

    #[error("Missing required SAM field: {0}")]
    MissingField(&'static str),

    #[error("Failed to parse integer in field '{field}': {source}")]
    ParseInt {
        field: &'static str,
        #[source]
        source: std::num::ParseIntError,
    },

    #[error("Invalid CIGAR character encountered: '{0}'")]
    InvalidCigar(char),

    #[error("Malformed SAM record: {0}")]
    MalformedRecord(String),
}

// ============================================================================
// Traits
// ============================================================================

/// Core SAM fields access
pub trait SamFields {
    fn qname(&self) -> &str;   // QNAME  (str): Query template NAME
    fn rname(&self) -> &str;   // RNAME  (str): Reference sequence NAME
    fn pos(&self)   -> u32;    // POS    (int): 1-based leftmost mapping POSition, 0 = unmapped
    fn mapq(&self)  -> u32;     // MAPQ   (int): MAPping Quality
    fn cigar(&self) -> &str;   // CIGAR  (str): CIGAR string
    fn rnext(&self) -> &str;   // RNEXT  (str): Reference name of the mate/next read
    fn pnext(&self) -> u32;    // PNEXT  (int): Position of the mate/next read
    fn tlen(&self)  -> i32;    // TLEN   (int): observed Template LENgth
    fn seq(&self)   -> &str;   // SEQ    (str): segment SEQuence
    fn qual(&self)  -> &str;   // QUAL   (str): ASCII of Phred-scaled base QUALity+33
}

/// Allow Flag based method implementation for SamStr and SamRecord
pub trait SamFlags: SamFields {
    fn flag(&self)             -> u16;
    fn is_proper(&self)        -> bool { (self.flag() & 0x0002) != 0 }
    fn is_mapped(&self)        -> bool { (self.flag() & 0x0004) == 0 }
    fn is_mate_mapped(&self)   -> bool { (self.flag() & 0x0008) == 0 }
    fn is_revcomp(&self)       -> bool { (self.flag() & 0x0010) != 0 }
    fn is_mate_revcomp(&self)  -> bool { (self.flag() & 0x0020) != 0 }
    fn is_read1(&self)         -> bool { (self.flag() & 0x0040) != 0 }
    fn is_read2(&self)         -> bool { (self.flag() & 0x0080) != 0 }
    fn is_secondary(&self)     -> bool { (self.flag() & 0x0100) != 0 }
    fn is_supplementary(&self) -> bool { (self.flag() & 0x0800) != 0 }
    fn is_primary(&self)       -> bool { !self.is_secondary() && !self.is_supplementary() }
}

/// Generic access to optional SAM tags
pub trait SamTags: SamFields {
    fn get_raw_tag(&self, tag: &str)   -> Option<&str>;
    fn get_int_tag(&self, tag: &str)   -> Option<i32> { self.get_raw_tag(tag).and_then(|val| val.parse::<i32>().ok()) }
    fn get_str_tag(&self, tag: &str)   -> Option<&str> { self.get_raw_tag(tag) }
    fn get_float_tag(&self, tag: &str) -> Option<f32> { self.get_raw_tag(tag).and_then(|val| val.parse::<f32>().ok()) }
}

/// Trait for parsing CIGAR string
pub trait CigarString: SamFields {
    /// Returns: Option<(AlignmentLength, TotalIndels, SoftClips, ReadLength)>
    fn parse_cigar(&self) -> Result<Option<(u32, u32, u32, u32)>, SamError>;
}

/// Blanket CigarString implementation for anything that implements SamFields.
impl <T: SamFields>CigarString for T {
    fn parse_cigar(&self) -> Result<Option<(u32, u32, u32, u32)>, SamError> {
        let cigar_str = self.cigar();
        if cigar_str == "*" || cigar_str.is_empty() {
            return Ok(None);
        }

        let mut align_len:   u32 = 0u32;
        let mut indels:      u32 = 0u32;
        let mut soft_clips:  u32 = 0u32;
        let mut read_len:    u32 = 0;
        let mut current_num: u32 = 0u32;

        for cigar_char in cigar_str.chars() {
            if cigar_char.is_ascii_digit() {
                current_num = current_num
                    .checked_mul(10)
                    .and_then(|n| n.checked_add(cigar_char.to_digit(10)?))
                    .unwrap_or(u32::MAX);
            } else {
                match cigar_char {
                    'M' | '=' | 'X' => {
                        align_len += current_num;
                        read_len += current_num;
                    }
                    'I' => {
                        align_len += current_num;
                        indels += current_num;
                        read_len += current_num;
                    }
                    'D' => {
                        align_len += current_num;
                        indels += current_num;
                    }
                    'S' => {
                        soft_clips += current_num;
                        read_len += current_num;
                    }
                    'N' | 'P' => {
                        // N is a reference skip (RNA-seq splicing)
                        // P is padding,
                    }
                    'H' => {
                        // Hard clips consume the original read, but are physically 
                        // removed from the SAM SEQ column.
                    }
                    _ => return Err(SamError::InvalidCigar(cigar_char)),
                }
                current_num = 0;
            }
        }
        
        Ok(Some((align_len, indels, soft_clips, read_len)))
    }
}

/// Trait for some custom alignment stats
pub trait AlignmentStats: CigarString + SamTags {
    fn calculate_as_al(&self) -> Result<Option<f32>, SamError>;                 // AS (alignment score) / AL (alignment length)
    fn calculate_alignment_length(&self) -> Result<Option<u32>, SamError>;      // from cigar string
    fn calculate_alignment_proportion(&self) -> Result<Option<f32>, SamError>;  // alignment length/read length
    fn calculate_alignment_accuracy(&self) -> Result<Option<f32>, SamError>;    // (alignment length - NM)/alignment length
}

/// Blanket AlignmentStats implementation for anything that implements SamFields, SamTags & CigarString.
impl <T: SamFields + SamTags + CigarString>AlignmentStats for T {
    
    /// Calculate AS/AL (Alignment Score divided by Alignment Length)
    fn calculate_as_al(&self) -> Result<Option<f32>, SamError> {
        let Some((align_len, _, _, _)) = self.parse_cigar()? else { return Ok(None) };
        let Some(as_score) = self.get_int_tag("AS") else { return Ok(None) };
        
        if align_len == 0 {
            Ok(None)
        } else {
            Ok(Some(as_score as f32 / align_len as f32))
        }
    }

    /// Get the alignment length from the CIGAR string
    fn calculate_alignment_length(&self) -> Result<Option<u32>, SamError> {
        Ok(self.parse_cigar()?.map(|(align_len, _, _, _)| align_len))
    }

    /// Calculate alignment proportion - % of read that is aligned
    fn calculate_alignment_proportion(&self) -> Result<Option<f32>, SamError> {
        let Some((align_len, _, _, read_len)) = self.parse_cigar()? else { return Ok(None) };
        if read_len == 0 {
            Ok(None)
        } else {
            Ok(Some(align_len as f32 / read_len as f32))
        }
    }

    /// Calculate alignment accuracy - i.e. % identity
    fn calculate_alignment_accuracy(&self) -> Result<Option<f32>, SamError> {
        let Some((align_len, _, _, _)) = self.parse_cigar()? else { return Ok(None) };
        let Some(nm) = self.get_int_tag("NM") else { return Ok(None) };
        
        if align_len == 0 {
            return Ok(None);
        }
        
        let matches = align_len.saturating_sub(nm as u32);
        Ok(Some((matches as f32 / align_len as f32) * 100.0))
    }

}

// ============================================================================
// Struct for Borrowed SAM record: SamStr (maybe change to BorrowedSam)
// ============================================================================


/// Struct representing a single SAM alignment record in a zero-allocation str format
pub struct SamStr<'a> {
    pub raw_line: &'a str,
}

impl<'a> SamStr<'a> {
    pub fn new(line: &'a str) -> Self {
        Self { raw_line: line }
    }
    
    pub fn to_fields(&self) -> impl Iterator<Item = &'a str> {
        self.raw_line.split('\t')
    }
    
    pub fn to_fastq_fields(&self) -> (&'a str, Cow<'a, str>, Cow<'a, str>) {
        let mut fields = self.to_fields();
        let qname: &str = fields.next().unwrap_or("*");
        let seq: &str   = fields.nth(8).unwrap_or("*");
        let qual: &str  = fields.next().unwrap_or("*");
        //(qname, seq, qual)
        if self.is_revcomp() { 
            let rc_seq: String = seq.reverse_complement();
            let rev_qual: String = qual.chars().rev().collect::<String>();
            (qname, Cow::Owned(rc_seq), Cow::Owned(rev_qual))
        } else {
            (qname, Cow::Borrowed(seq), Cow::Borrowed(qual))
        }
    }

    pub fn seq_length(&self) -> usize {
        self.seq().len()
    }
}

impl<'a> SamFields for SamStr<'a> {
    fn qname(&self) -> &str { self.to_fields().next().unwrap_or("*") }
    fn rname(&self) -> &str { self.to_fields().nth(2).unwrap_or("*") }
    fn pos(&self)   -> u32 { self.to_fields().nth(3).and_then(|s| s.parse().ok()).unwrap_or(0) }
    fn mapq(&self)  -> u32 { self.to_fields().nth(4).and_then(|s| s.parse().ok()).unwrap_or(0) }
    fn cigar(&self) -> &str { self.to_fields().nth(5).unwrap_or("*") }
    fn rnext(&self) -> &str { self.to_fields().nth(6).unwrap_or("*") }
    fn pnext(&self) -> u32 { self.to_fields().nth(7).and_then(|s| s.parse().ok()).unwrap_or(0) }
    fn tlen(&self)  -> i32 { self.to_fields().nth(8).and_then(|s| s.parse().ok()).unwrap_or(0) }
    fn seq(&self)   -> &str { self.to_fields().nth(9).unwrap_or("*") }
    fn qual(&self)  -> &str { self.to_fields().nth(10).unwrap_or("*") }
}

impl<'a> SamFlags for SamStr<'a> {
    fn flag(&self) -> u16 {
        self.to_fields().nth(1).and_then(|s| s.parse().ok()).unwrap_or(0)
    }
}

impl<'a> SamTags for SamStr<'a> {
    fn get_raw_tag(&self, tag: &str) -> Option<&str> {
        if tag.len() != 2 { return None; }
        self.to_fields().skip(11).find_map(|field| {
            let bytes: &[u8] = field.as_bytes();
            if field.starts_with(tag) && bytes.get(2) == Some(&b':') && bytes.get(4) == Some(&b':') {
                Some(&field[5..])
            } else {
                None
            }
        })
    }
}

// ============================================================================
// Struct for Owned SAM record: SamRecord (maybe change to OwnedSam)
// ============================================================================

/// Struct representing a fully parsed, owned SAM alignment record.
#[derive(Debug, Clone)]
pub struct SamRecord {
    pub qname: String,
    pub flag:  u16,
    pub rname: String,
    pub pos:   u32,
    pub mapq:  u32,
    pub cigar: String,
    pub rnext: String,
    pub pnext: u32,
    pub tlen:  i32,
    pub seq:   String,
    pub qual:  String,
    pub tags:  Vec<String>,
}

impl SamRecord {
    pub fn from_str(line: &str) -> Result<Self, SamError> {
        SamRecord::try_from(&SamStr::new(line))
    }
    
    pub fn seq_length(&self) -> usize {
        self.seq.len()
    }
    
    pub fn passes_mapq(&self, min_mapq: u32) -> bool {
        self.mapq >= min_mapq
    }
    
    pub fn into_fastq<T: FastqFormatter>(self) -> FastqRecord<T> {
        let (final_seq, final_qual) = if self.is_revcomp() {
            let rc_seq: String = self.seq.reverse_complement();
            let rev_qual: String = self.qual.chars().rev().collect::<String>();
            (rc_seq, rev_qual)
        } else {
            (self.seq, self.qual)
        };

        FastqRecord::new(self.qname, final_seq, final_qual)
    }
}

impl<'a> TryFrom<&SamStr<'a>> for SamRecord {
    type Error = SamError;

    fn try_from(sam: &SamStr<'a>) -> Result<Self, Self::Error> {
        let mut fields = sam.to_fields();

        let qname: String     = fields.next().ok_or(SamError::MissingField("QNAME"))?.to_string();
        let flag: u16         = fields.next().and_then(|s| s.parse().ok()).unwrap_or(0);
        let rname: String     = fields.next().unwrap_or("*").to_string();
        let pos: u32          = fields.next().and_then(|s| s.parse().ok()).unwrap_or(0);
        let mapq: u32         = fields.next().and_then(|s| s.parse().ok()).unwrap_or(0);
        let cigar: String     = fields.next().unwrap_or("*").to_string();
        let rnext: String     = fields.next().unwrap_or("*").to_string();
        let pnext: u32        = fields.next().and_then(|s| s.parse().ok()).unwrap_or(0);
        let tlen: i32         = fields.next().and_then(|s| s.parse().ok()).unwrap_or(0);
        let seq: String       = fields.next().unwrap_or("*").to_string();
        let qual: String      = fields.next().unwrap_or("*").to_string();
        let tags: Vec<String> = fields.map(|s| s.to_string()).collect();

        Ok(SamRecord {
            qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, tags
        })
    }
}

impl SamFields for SamRecord {
    fn qname(&self) -> &str { &self.qname }
    fn rname(&self) -> &str { &self.rname }
    fn pos(&self)   -> u32  { self.pos }
    fn mapq(&self)  -> u32  { self.mapq }
    fn cigar(&self) -> &str { &self.cigar }
    fn rnext(&self) -> &str { &self.rnext }
    fn pnext(&self) -> u32  { self.pnext }
    fn tlen(&self)  -> i32  { self.tlen }
    fn seq(&self)   -> &str { &self.seq }
    fn qual(&self)  -> &str { &self.qual }
}

impl SamFlags for SamRecord {
    fn flag(&self) -> u16 { self.flag }
}

impl SamTags for SamRecord {
    fn get_raw_tag(&self, tag: &str) -> Option<&str> {
        if tag.len() != 2 { return None; }
        self.tags.iter().find_map(|field| {
            let bytes: &[u8] = field.as_bytes();
            if field.starts_with(tag) && bytes.get(2) == Some(&b':') && bytes.get(4) == Some(&b':') {
                Some(&field[5..])
            } else {
                None
            }
        })
    }
}

// ============================================================================
// Struct for SAM Reader
// ============================================================================

/// The main reader for SAM files.
pub struct SamReader {
    pub reader:  ReaderType,
    pub headers: Vec<String>,
}

impl SamReader {
    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self, SamError> {
        let file = File::open(path)?;
        Ok(SamReader {
            reader: ReaderType::File(BufReader::new(file)),
            headers: Vec::new(),
        })
    }

    pub fn from_stdin() -> Self {
        SamReader {
            reader:  ReaderType::Stdin(BufReader::new(io::stdin())),
            headers: Vec::new(),
        }
    }
}

impl Iterator for SamReader {
    type Item = Result<SamRecord, SamError>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut line = String::new();
        loop {
            line.clear();
            match self.reader.read_line(&mut line) {
                Ok(0) => return None,
                Ok(_) => {
                    let trimmed = line.trim_end();
                    if trimmed.is_empty() { continue; }
                    if trimmed.starts_with('@') {
                        self.headers.push(trimmed.to_string());
                        continue;
                    }
                    return Some(SamRecord::from_str(trimmed));
                }
                Err(e) => return Some(Err(SamError::Io(e))),
            }
        }
    }
}

/// A reader optimized for parallel processing chunk by chunk via Rayon.
pub struct ChunkySamReader<R: BufRead> {
    reader:         R,
    pub headers:    Vec<String>,
    pub chunk_size: usize,
}

impl<R: BufRead> ChunkySamReader<R> {
    pub fn new(reader: R, chunk_size: usize) -> Self {
        Self {
            reader,
            headers: Vec::new(),
            chunk_size,
        }
    }

    pub fn next_chunk(&mut self) -> Result<Option<Vec<String>>, SamError> {
        let mut chunk = Vec::with_capacity(self.chunk_size);
        let mut line_buffer = String::with_capacity(1024);

        while chunk.len() < self.chunk_size {
            line_buffer.clear();
            match self.reader.read_line(&mut line_buffer) {
                Ok(0) => break, // EOF
                Ok(_) => {
                    let trimmed = line_buffer.trim_end();
                    if trimmed.is_empty() { continue; }

                    if trimmed.starts_with('@') {
                        self.headers.push(trimmed.to_string());
                    } else {
                        chunk.push(trimmed.to_string());
                    }
                }
                Err(e) => return Err(SamError::Io(e)),
            }
        }

        if chunk.is_empty() {
            Ok(None)
        } else {
            Ok(Some(chunk))
        }
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    const TEST_SAM_LINE: &str = "read_1\t99\tchr1\t10000\t60\t10M2I3S\t=\t10200\t210\tACGTACGTACIIISS\tIIIIIIIIIIIIIII\tNM:i:2\tAS:i:10\tXZ:Z:dummy";

    #[test]
    fn test_sam_str_parsing() {
        let sam: SamStr<'_> = SamStr::new(TEST_SAM_LINE);
        assert_eq!(sam.qname(), "read_1");
        assert_eq!(sam.flag(), 99);
        assert_eq!(sam.rname(), "chr1");
        assert_eq!(sam.pos(), 10000);
        assert_eq!(sam.mapq(), 60);
        assert_eq!(sam.cigar(), "10M2I3S");
        assert_eq!(sam.seq(), "ACGTACGTACIIISS");
        assert_eq!(sam.qual(), "IIIIIIIIIIIIIII");
    }

    #[test]
    fn test_sam_tags() {
        let sam: SamStr<'_> = SamStr::new(TEST_SAM_LINE);
        assert_eq!(sam.get_int_tag("NM"), Some(2));
        assert_eq!(sam.get_int_tag("AS"), Some(10));
        assert_eq!(sam.get_str_tag("XZ"), Some("dummy"));
        assert_eq!(sam.get_str_tag("XX"), None);
    }

    #[test]
    fn test_cigar_parsing() {
        let sam: SamStr<'_> = SamStr::new(TEST_SAM_LINE);
        // 10M + 2I = 12 align_len. Indels = 2. SoftClips = 3. Read_len = 10 + 2 + 3 = 15.
        let parsed: Option<(u32, u32, u32, u32)> = sam.parse_cigar().unwrap();
        assert_eq!(parsed, Some((12u32, 2u32, 3u32, 15u32)));
    }

    #[test]
    fn test_alignment_stats() {
        let sam: SamStr<'_> = SamStr::new(TEST_SAM_LINE);
        // AS/AL -> 10 / 12 = 0.8333...
        let as_al: Option<f32> = sam.calculate_as_al().unwrap();
        assert!((as_al.unwrap() - 0.8333).abs() < 0.001);

        // Accuracy -> (12 - 2) / 12 = 83.33%
        let acc: Option<f32> = sam.calculate_alignment_accuracy().unwrap();
        assert!((acc.unwrap() - 83.33).abs() < 0.01);
    }

    #[test]
    fn test_sam_record_conversion() {
        let record = SamRecord::from_str(TEST_SAM_LINE).unwrap();
        assert_eq!(record.qname, "read_1");
        assert!(record.is_read1());
        assert!(record.is_proper());
        assert!(!record.is_revcomp());
        assert!(record.is_mate_revcomp());
    }
}