//! n2core/src/sequence.rs

//use std::fmt;

/// A validated Nucleotide sequence 
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NucleotideSequence(pub Vec<u8>);

/// A validated DNA sequence 
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DnaSequence(pub Vec<u8>);

/// A validated RNA sequence 
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RnaSequence(pub Vec<u8>);

/// A validated sequence of amino acids.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProteinSequence(pub Vec<u8>);

/// Represents a 3-base sequence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Codon(pub [u8; 3]);

/// Represents a single amino acid character.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AminoAcid(pub u8);


pub fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'a' => 't',
            'C' => 'G',
            'c' => 'g',
            'G' => 'C',
            'g' => 'c',
            'T' => 'A',
            't' => 'a',
            'N' => 'N',
            'n' => 'n',
            _ => c, // Pass through any unexpected characters
        })
        .collect()
}