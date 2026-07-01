//! n2bio/pfqsim/src/genome.rs
//! 

use std::io::{ self, Error, ErrorKind };
use rand::RngExt; 
use rand::distr::{ Distribution, weighted::WeightedIndex };

use n2core::fasta::{ FastaReader, FastaRecord };
use n2core::readers::ReaderType;

// ============================================================================
// Contig
// ============================================================================

pub(crate) struct Contig {
    pub(crate) id: String,
    pub(crate) seq: Vec<u8>, 
}

// ============================================================================
// Reference genome
// ============================================================================

pub(crate) struct ReferenceGenome {
    pub(crate) contigs: Vec<Contig>,
    pub(crate) selector: WeightedIndex<usize>,
}

impl ReferenceGenome {
    /// Inspects a FASTA file path to calculate total length (without N's) and contig counts 
    pub(crate) fn parse_fasta_metrics(path: &str) -> io::Result<(usize, usize)> {
        let reader: FastaReader<ReaderType> = FastaReader::open(path)?;
        
        let mut clean_length: usize = 0;
        let mut contig_count: usize = 0;

        for record_result in reader {
            let record: FastaRecord = record_result.map_err(|e| {
                io::Error::new(io::ErrorKind::InvalidData, format!("FASTA parsing error: {}", e))
            })?;
            
            contig_count += 1;

            for &byte in record.seq.as_bytes() {
                if byte != b'N' && byte != b'n' {
                    clean_length += 1;
                }
            }
        }

        Ok((clean_length, contig_count))
    }

    /// Loads the FASTA, handles circularity, splits at 'N's, and builds weights
    pub(crate) fn load(
        reader: FastaReader<ReaderType>, 
        min_length: usize, 
        is_circular: bool,
    ) -> io::Result<Self> {
        let mut contigs: Vec<Contig> = Vec::new();
        let mut weights: Vec<usize>  = Vec::new();
        let mut total_bases: usize   = 0;

        for record_result in reader {
            let record: FastaRecord  = record_result?;
            let mut raw_seq: Vec<u8> = record.seq.into_bytes();

            // 1. Handle Circularity before splitting
            if is_circular && raw_seq.len() >= min_length {
                let bridge: Vec<u8> = raw_seq[0..min_length].to_vec();
                raw_seq.extend(bridge);
            }

            // 2. Split by 'N' or 'n' to remove from the sampling pool
            for chunk in raw_seq.split(|&b| b == b'N' || b == b'n') {
                let seq_len: usize = chunk.len();

                if seq_len >= min_length {
                    total_bases += seq_len;
                    contigs.push(Contig {
                        id: record.id.clone(), 
                        seq: chunk.to_vec(),
                    });
                    weights.push(seq_len);
                }
            }
        }

        if contigs.is_empty() {
            return Err(Error::new(
                ErrorKind::InvalidData,
                format!("No sequences (or 'N'-free chunks) are >= {} bp.", min_length),
            ));
        }

        println!("Loaded {} valid contig chunks ({} clean bases total)", contigs.len(), total_bases);

       
        let selector: WeightedIndex<usize> = WeightedIndex::new(weights)
            .map_err(|_| Error::new(ErrorKind::Other, "Failed to build weighted index"))?;

        Ok(Self { contigs, selector })
    }

    /// Randomly selects a contig (weighted by length) and returns the ID and a raw byte slice
    pub(crate) fn sample_slice<'a, R: rand::Rng + ?Sized>(
        &'a self,
        rng: &mut R,
        insert_size: usize,
        buffer_size: usize,
    ) -> (&'a str, &'a [u8]) {
        // Sample a slice from distribution using rng
        let contig_idx: usize   = self.selector.sample(rng);
        let contig: &Contig     = &self.contigs[contig_idx];
        let slice_length: usize = insert_size + buffer_size;
        let max_start: usize    = contig.seq.len().saturating_sub(slice_length);
        let start_pos: usize    = if max_start == 0 { 0 } else { rng.random_range(0..=max_start) };
        let slice_end: usize    = std::cmp::min(start_pos + slice_length, contig.seq.len());

        (&contig.id, &contig.seq[start_pos..slice_end])
    }
}