//! n2core/src/sequence.rs
//! Core traits and utilities for DNA sequence manipulation.

use thiserror::Error;

// ============================================================================
// Errors
// ============================================================================

/// Custom errors for sequence processing
#[derive(Error, Debug, PartialEq)]
pub enum SequenceError {
    #[error("K-mer size cannot be zero.")]
    ZeroLengthKmer,

    #[error("DNA string contains non-ASCII characters.")]
    NonAsciiSequence,
}

// ============================================================================
// Traits
// ============================================================================

/// Core functions for DNA sequences
pub trait DnaSequence {
    /// Owned type (String or Vec<u8>)
    type OwnedSeq;
    
    /// Borrowed slice (&str or &[u8])
    type SeqSlice<'a> where Self: 'a;

    /// Returns the reverse complement of the sequence
    fn reverse_complement(&self) -> Self::OwnedSeq;
    
    /// Returns an iterator yielding consecutive k-mers of length `k`
    /// 
    /// # Errors
    /// Returns `SequenceError::ZeroLengthKmer` if `k` is 0.
    /// Returns `SequenceError::NonAsciiSequence` if implemented on `str` and the text is not valid ASCII.
    fn to_kmers(&self, k: usize) -> Result<impl Iterator<Item = Self::SeqSlice<'_>>, SequenceError>;
}

// ============================================================================
// Implementation for str
// ============================================================================

impl DnaSequence for str {
    type OwnedSeq = String;
    type SeqSlice<'a> = &'a str;

    fn reverse_complement(&self) -> Self::OwnedSeq {
        self.chars()
            .rev()
            .map(|c| match c {
                'A' => 'T', 'a' => 't',
                'C' => 'G', 'c' => 'g',
                'G' => 'C', 'g' => 'c',
                'T' => 'A', 't' => 'a',
                'N' => 'N', 'n' => 'n',
                _ => c, 
            })
            .collect()
    }

    fn to_kmers(&self, k: usize) -> Result<impl Iterator<Item = Self::SeqSlice<'_>>, SequenceError> {
        if k == 0 {
            return Err(SequenceError::ZeroLengthKmer);
        }
        
        if !self.is_ascii() {
            return Err(SequenceError::NonAsciiSequence);
        }

        Ok(self.as_bytes().windows(k).map(|w| {
            unsafe { std::str::from_utf8_unchecked(w) }
        }))
    }
}

// ============================================================================
// Implementation for u8
// ============================================================================

impl DnaSequence for [u8] {
    type OwnedSeq = Vec<u8>;
    type SeqSlice<'a> = &'a [u8];

    fn reverse_complement(&self) -> Self::OwnedSeq {
        self.iter()
            .rev()
            .map(|&c| match c {
                b'A' => b'T', b'a' => b't',
                b'C' => b'G', b'c' => b'g',
                b'G' => b'C', b'g' => b'c',
                b'T' => b'A', b't' => b'a',
                b'N' => b'N', b'n' => b'n',
                _ => c, // Keep other bytes as-is
            })
            .collect()
    }

    fn to_kmers(&self, k: usize) -> Result<impl Iterator<Item = Self::SeqSlice<'_>>, SequenceError> {
        if k == 0 {
            return Err(SequenceError::ZeroLengthKmer);
        }
        
        Ok(self.windows(k))
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement_str() {
        let seq: &str = "ACGTacgtNn-";
        let rc: String = seq.reverse_complement();
        assert_eq!(rc, "-nNacgtACGT");
    }

    #[test]
    fn test_reverse_complement_u8() {
        let seq: &[u8; 11] = b"ACGTacgtNn-";
        let rc: Vec<u8> = seq.reverse_complement();
        assert_eq!(rc, b"-nNacgtACGT");
    }

    #[test]
    fn test_to_kmers_str() {
        let seq: &str = "ACGT";
        let mut kmers = seq.to_kmers(3).unwrap();
        assert_eq!(kmers.next(), Some("ACG"));
        assert_eq!(kmers.next(), Some("CGT"));
        assert_eq!(kmers.next(), None);
    }

    #[test]
    fn test_to_kmers_u8() {
        let seq: &[u8; 4] = b"ACGT";
        let mut kmers = seq.to_kmers(3).unwrap();
        assert_eq!(kmers.next(), Some(&b"ACG"[..]));
        assert_eq!(kmers.next(), Some(&b"CGT"[..]));
        assert_eq!(kmers.next(), None);
    }

    #[test]
    fn test_to_kmers_zero_length() {
        let seq_str: &str = "ACGT";
        assert!(matches!(
            seq_str.to_kmers(0),
            Err(SequenceError::ZeroLengthKmer)
        ));

        let seq_u8: &[u8; 4] = b"ACGT";
        assert!(matches!(
            seq_u8.to_kmers(0),
            Err(SequenceError::ZeroLengthKmer)
        ));
    }

    #[test]
    fn test_to_kmers_non_ascii_str() {
        // Includes a multi-byte character (🧬)
        let seq: &str = "ACG🧬T";
        assert!(matches!(
            seq.to_kmers(2),
            Err(SequenceError::NonAsciiSequence)
        ));
    }
}