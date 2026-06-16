//! n2bio/pfqsim/src/mutate.rs
//! 

use rand::RngExt;
use rand::rngs::SmallRng;

// ============================================================================
// Data Structures
// ============================================================================

pub struct Mutator {
    pub sub_rate: f64,
    pub indel_rate: f64,
}

pub struct MutationStats {
    pub sequence: Vec<u8>,
    pub subs: usize,
    pub insertions: usize,
    pub deletions: usize,
}

// ============================================================================
// Implementation
// ============================================================================

impl Mutator {
    pub fn new(sub_rate: f64, indel_rate: f64) -> Self {
        Self {
            sub_rate,
            indel_rate,
        }
    }

    /// Mutates a raw reference slice in a single pass
    pub fn mutate(
        &self,
        ref_slice: &[u8],
        target_length: usize,
        rng: &mut SmallRng,
    ) -> MutationStats {
        let mut sequence: Vec<u8> = Vec::with_capacity(target_length);
        let mut subs: usize = 0;
        let mut insertions: usize = 0;
        let mut deletions: usize = 0;
        let mut ref_idx: usize = 0;

        while sequence.len() < target_length && ref_idx < ref_slice.len() {
            let current_base: u8 = ref_slice[ref_idx];

            // 1. Roll for Indels
            if self.indel_rate > 0.0 && rng.random_bool(self.indel_rate) {
                if rng.random_bool(0.5) {
                    // Insertion: add a random base, do not advance the reference index
                    sequence.push(random_base(rng));
                    insertions += 1;
                    continue; 
                } else {
                    // Deletion: advance the reference index, do not add to sequence
                    ref_idx += 1;
                    deletions += 1;
                    continue;
                }
            }

            // 2. Roll for Substitutions
            if self.sub_rate > 0.0 && rng.random_bool(self.sub_rate) {
                sequence.push(random_mutated_base(rng, current_base));
                subs += 1;
            } else {
                // Force soft-masked lowercase reference sequence to uppercase
                sequence.push(current_base.to_ascii_uppercase());
            }

            ref_idx += 1;
        }

        // Fallback: Pad with 'N's if severe deletions exhausted the buffer
        //while sequence.len() < target_length {
        //    sequence.push(b'N');
        
        // Fallback: pad with random nucleotide
        while sequence.len() < target_length {
            sequence.push(random_base(rng));
        }

        MutationStats {
            sequence,
            subs,
            insertions,
            deletions,
        }
    }
}

// ============================================================================
// Helper Functions
// ============================================================================

const BASES: &[u8] = b"ACGT";

#[inline]
fn random_base(rng: &mut SmallRng) -> u8 {
    BASES[rng.random_range(0..4)]
}

#[inline]
fn random_mutated_base(rng: &mut SmallRng, original: u8) -> u8 {
    loop {
        let new_base = random_base(rng);
        // Compare case in case FASTA was lowercase
        if new_base.to_ascii_uppercase() != original.to_ascii_uppercase() {
            return new_base;
        }
    }
}