//! n2core/src/kmer.rs

use std::hash::{Hash, Hasher};
use std::collections::hash_map::DefaultHasher;
use crate::sequence::DnaSequence;

pub trait KmerHash {
    /// Returns the canonical version of the k-mer for strand-agnostic graph building
    fn canonical(&self) -> Vec<u8>;
    
    /// Generate a standard u64 hash of the k-mer
    fn kmer_hash(&self) -> u64;
    
    /// Helper to determine if the k-mer is a minimizer (e.g., lower hash value)
    fn is_minimizer_against(&self, other: &Self) -> bool;
}

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

pub trait KmerEncoding {
    /// Compresse a kmer (byte slice) into a u64 integer
    fn encode_to_u64(&self) -> u64;
    
    /// Decode a u64 integer back into a DNA sequence of length `k`
    fn decode_from_u64(encoded: u64, k: usize) -> Vec<u8>;

    /// Calculate the canonical (lower value) encoding
    fn canonical_u64(&self) -> u64;
}

impl KmerEncoding for [u8] {
    fn encode_to_u64(&self) -> u64 {
        assert!(self.len() <= 32, "k-mer must be 32 bases or fewer to fit in a u64");
        
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
        encoded
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

    fn canonical_u64(&self) -> u64 {
        let fwd_u64 = self.encode_to_u64();
        let rc_u64 = self.reverse_complement().encode_to_u64();
        
        // Return the lesser one
        std::cmp::min(fwd_u64, rc_u64)
    }
}

