//! n2core/src/kmer.rs
//! Utilities for k-mer hashing, encoding, and strand orientation.

use std::collections::hash_map::DefaultHasher;
use std::collections::HashSet;
use std::hash::{Hash, Hasher};
use thiserror::Error;

use crate::sequence::{DnaSequence, SequenceError};

// ============================================================================
// Errors
// ============================================================================

#[derive(Error, Debug, PartialEq)]
pub enum KmerError {
    #[error("K-mer length {0} exceeds the maximum allowed (32) for u64 encoding.")]
    KmerTooLong(usize),

    #[error("Sequence error: {0}.")]
    Sequence(#[from] SequenceError),
}

// ============================================================================
// Traits
// ============================================================================

/// Hash based kmer encoding
pub trait KmerHash {
    /// Returns the canonical version of the k-mer
    fn canonical(&self) -> Vec<u8>;
    
    /// Generate a standard u64 hash of the k-mer
    fn kmer_hash(&self) -> u64;
    
    /// Helper to determine if the k-mer is a minimizer (e.g., lower hash value)
    fn is_minimizer_against(&self, other: &Self) -> bool;
}

/// Two-bit kmer encoding
pub trait KmerEncoding {
    /// Compresse a kmer (byte slice) into a u64 integer
    fn encode_to_u64(&self) -> Result<u64, KmerError>;
    
    /// Decode a u64 integer back into a DNA sequence of length `k`
    fn decode_from_u64(encoded: u64, k: usize) -> Vec<u8>;

    /// Calculate the canonical (lower value) encoding
    fn canonical_u64(&self) -> Result<u64, KmerError>;
}

// ============================================================================
// Implementations
// ============================================================================

impl KmerHash for [u8] {
    fn canonical(&self) -> Vec<u8> {
        let rev_comp: Vec<u8> = self.reverse_complement();
        // Return whichever sequence is lexicographically first
        if self < &rev_comp[..] {
            self.to_vec()
        } else {
            rev_comp
        }
    }

    fn kmer_hash(&self) -> u64 {
        let mut hasher: DefaultHasher = DefaultHasher::new();
        self.canonical().hash(&mut hasher);
        hasher.finish()
    }

    fn is_minimizer_against(&self, other: &Self) -> bool {
        self.kmer_hash() < other.kmer_hash()
    }
}

impl KmerEncoding for [u8] {
    fn encode_to_u64(&self) -> Result<u64, KmerError> {
        if self.len() > 32 {
            return Err(KmerError::KmerTooLong(self.len()));
        }
        let mut encoded: u64 = 0u64;
        for &base in self {
            encoded <<= 2; // Shift left by 2 bits for the next base
            encoded |= match base {
                b'A' | b'a' => 0b00,
                b'C' | b'c' => 0b01,
                b'G' | b'g' => 0b10,
                b'T' | b't' => 0b11,
                _           => 0b00, // 'N' is coerced to 'A' in strict 2-bit encoding
            };
        }
        Ok(encoded)
    }

    fn decode_from_u64(mut encoded: u64, k: usize) -> Vec<u8> {
        let mut sequence: Vec<u8> = vec![0u8; k];
        
        // Decode backwards (from the last base to the first)
        for i in (0..k).rev() {
            let base_bits: u64 = encoded & 0b11; // Mask all but the last 2 bits
            sequence[i] = match base_bits {
                0b00 => b'A',
                0b01 => b'C',
                0b10 => b'G',
                0b11 => b'T',
                _    => unreachable!(),
            };
            encoded >>= 2; // Shift right by 2 bits for the next iteration
        }
        
        sequence
    }

    fn canonical_u64(&self) -> Result<u64, KmerError> {
        let fwd_u64: u64 = self.encode_to_u64()?;
        let rc_u64: u64 = self.reverse_complement().encode_to_u64()?;
        
        // Return the lesser one
        Ok(std::cmp::min(fwd_u64, rc_u64))
    }
}


// ============================================================================
// Kmer based strand orientor
// ============================================================================

pub struct StrandOrientor {
    // A set of strictly directional (forward) 2-bit encoded k-mers from the reference
    reference_kmers: HashSet<u64>,
    k: usize,
}

impl StrandOrientor {
    // Initializes the orientor using the reference backbone
    pub fn new(reference: &[u8], k: usize) -> Result<Self, KmerError> {
        let mut reference_kmers: HashSet<u64> = HashSet::new();
        let kmers = reference.to_kmers(k)?;
        for kmer in kmers {
            reference_kmers.insert(kmer.encode_to_u64()?);
        }

        Ok(Self { reference_kmers, k })
    }

    /// Test a sequence and return it correctly oriented to the reference strand
    pub fn orient(&self, sequence: &[u8]) -> Result<Vec<u8>, KmerError> {
        let mut fwd_hits: i32 = 0;
        let mut rc_hits: i32 = 0;

        let kmers = sequence.to_kmers(self.k)?;

        // Single loop: test both forward and RC encoded k-mers against the reference
        for kmer in kmers {
            let fwd_encoded: u64 = kmer.encode_to_u64()?;
            if self.reference_kmers.contains(&fwd_encoded) {
                fwd_hits += 1;
            }

            let rc_encoded: u64 = kmer.reverse_complement().encode_to_u64()?;
            if self.reference_kmers.contains(&rc_encoded) {
                rc_hits += 1;
            }
        }

        // Return whichever sequence mapped better to the reference strand.
        if rc_hits > fwd_hits {
            Ok(sequence.reverse_complement())
        } else {
            Ok(sequence.to_vec())
        }
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_decode_u64() {
        let seq: &[u8; 8] = b"ACGTACGT";
        let encoded: u64 = seq.encode_to_u64().unwrap();
        let decoded: Vec<u8> = <[u8]>::decode_from_u64(encoded, seq.len());
        assert_eq!(seq.as_ref(), decoded.as_slice());
    }

    #[test]
    fn test_encode_too_long() {
        // 33 bases (max is 32)
        let seq: &[u8; 33] = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"; 
        let err: KmerError = seq.encode_to_u64().unwrap_err();
        assert_eq!(err, KmerError::KmerTooLong(33));
    }

    #[test]
    fn test_canonical_u64() {
        let seq1: &[u8; 4] = b"ACGT";
        let seq2: Vec<u8> = b"ACGT".reverse_complement();
        let can1: u64 = seq1.canonical_u64().unwrap();
        let can2: u64 = seq2.canonical_u64().unwrap();
        assert_eq!(can1, can2);
    }

    #[test]
    fn test_strand_orientor() {
        // Test 1: Sequence is already in the forward orientation
        let reference: &[u8; 12] = b"ATGCATGCATGC";
        let orientor: StrandOrientor = StrandOrientor::new(reference, 4).unwrap();
        let fwd_query: &[u8; 6] = b"ATGCAT";
        let oriented_fwd: Vec<u8> = orientor.orient(fwd_query).unwrap();
        assert_eq!(oriented_fwd, b"ATGCAT");

        // Test 2: Sequence is reverse complemented and needs to be flipped
        let pure_fwd: &[u8; 7] = b"GATTACA"; 
        let pure_rc: Vec<u8> = pure_fwd.reverse_complement(); // b"TGTAATC"
        let orientor2: StrandOrientor = StrandOrientor::new(pure_fwd, 4).unwrap();
        let result = orientor2.orient(&pure_rc).unwrap();
        assert_eq!(result, pure_fwd);
    }
}