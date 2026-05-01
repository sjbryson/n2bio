//! n2core/src/sequence.rs

use std::hash::{Hash, Hasher};
use std::collections::hash_map::DefaultHasher;


pub trait DnaSequence {
    /// Owned type returned by reverse complement (String or Vec<u8>)
    type OwnedSeq;
    
    /// Borrowed slice yielded by the kmer iterator (&str or &[u8])
    type SeqSlice<'a> where Self: 'a;

    /// Returns the reverse complement of the sequence
    fn reverse_complement(&self) -> Self::OwnedSeq;
    
    /// Returns an iterator yielding consecutive k-mers of length `k`
    fn to_kmers(&self, k: usize) -> impl Iterator<Item = Self::SeqSlice<'_>>;
}

// ---------------------------------------------------------
// Implementation for OwnedSeq
// ---------------------------------------------------------
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

    fn to_kmers(&self, k: usize) -> impl Iterator<Item = Self::SeqSlice<'_>> {
        // Rust strings don't have a `.windows()` method natively because UTF-8 
        // characters can be variable length. Since DNA is strictly ASCII, 
        // we can safely drop down to bytes, window them, and convert back to &str.
        self.as_bytes()
            .windows(k)
            .map(|w| std::str::from_utf8(w).expect("Invalid UTF-8 in DNA string"))
    }
}

// ---------------------------------------------------------
// Implementation for SeqSlice
// ---------------------------------------------------------
impl DnaSequence for [u8] {
    type OwnedSeq = Vec<u8>;
    type SeqSlice<'a> = &'a [u8];

    fn reverse_complement(&self) -> Self::OwnedSeq {
        self.iter()
            .rev()
            .map(|&c| match c {
                b'A' | b'a' => b'T',
                b'T' | b't' => b'A',
                b'C' | b'c' => b'G',
                b'G' | b'g' => b'C',
                b'N' | b'n' => b'N',
                _ => c, // Keep other bytes as-is
            })
            .collect()
    }

    fn to_kmers(&self, k: usize) -> impl Iterator<Item = Self::SeqSlice<'_>> {
        self.windows(k)
    }
}


pub trait Kmer {
    /// Returns the canonical version of the k-mer for strand-agnostic graph building
    fn canonical(&self) -> Vec<u8>;
    
    /// Generates a standard u64 hash of the k-mer
    fn kmer_hash(&self) -> u64;
    
    /// Helper to determine if this k-mer is a minimizer (e.g., lower hash value)
    fn is_minimizer_against(&self, other: &Self) -> bool;
}

impl Kmer for [u8] {
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
        // In practice, for bioinformatics, you might want to swap DefaultHasher 
        // with something faster like ahash or wyhash later.
        self.canonical().hash(&mut hasher);
        hasher.finish()
    }

    fn is_minimizer_against(&self, other: &Self) -> bool {
        self.kmer_hash() < other.kmer_hash()
    }
}



