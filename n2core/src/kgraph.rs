//! n2core/src/kgraph.rs

use std::fs::File;
use std::io::{self, BufReader, BufWriter};
use serde::{Deserialize, Serialize};
use petgraph::graph::{DiGraph, NodeIndex};
use std::collections::HashMap;
use std::collections::HashSet;
use crate::sequence::DnaSequence;
use crate::kmer::KmerEncoding;
use crate::fasta::FastaReader;


/// Data stored in each node
#[derive(Serialize, Deserialize)]
pub struct KmerNode {
    pub canonical_u64: u64,
    pub frequency: u32,
    // Use HashSet<String> to track genomes associated with each kmer?
}
#[derive(Serialize, Deserialize)]
pub struct PanGenomeGraph {
    /// Petgraph storing nodes and edges
    pub graph: DiGraph<KmerNode, u32>,
    /// Lookup index mapping canonical k-mers to nodes
    pub node_map: HashMap<u64, NodeIndex>,
    /// K-mer size used for the graph
    pub k: usize,
}

impl PanGenomeGraph {
    pub fn new(k: usize) -> Self {
        Self {
            graph: DiGraph::new(),
            node_map: HashMap::new(),
            k,
        }
    }

    /// Incorporates a new sequence in the graph
    pub fn add_sequence(&mut self, seq: &str) {
        let mut prev_node: Option<NodeIndex> = None;

        // Calling .as_bytes(), uses the &[u8] implementation of DnaSequence
        for kmer in seq.as_bytes().to_kmers(self.k) {
            let canonical: u64 = kmer.canonical_u64();

            // 1. Get the existing node, or create a new one
            let current_node: NodeIndex = if let Some(&idx) = self.node_map.get(&canonical) {
                self.graph[idx].frequency += 1;
                idx
            } else {
                let new_node: KmerNode = KmerNode {
                    canonical_u64: canonical.clone(),
                    frequency: 1,
                };
                let idx: NodeIndex = self.graph.add_node(new_node);
                self.node_map.insert(canonical, idx);
                idx
            };

            // 2. Connect to the previous node (if this isn't the first k-mer)
            if let Some(prev) = prev_node {
                // If the edge exists, increment its frequency. Otherwise, create it.
                if let Some(edge_idx) = self.graph.find_edge(prev, current_node) {
                    self.graph[edge_idx] += 1;
                } else {
                    self.graph.add_edge(prev, current_node, 1);
                }
            }

            // Move the window forward
            prev_node = Some(current_node);
        }
    }

    /// Serialize the graph to a binary file
    /// Example: Build and save
    /// let graph = PanGenomeGraph::from_fastas("ref.fa", "assemblies.fa", 31)?;
    /// graph.save_to_file("virus_pangenome.bin")?;
    /// 
    pub fn save_to_file(&self, path: &str) -> io::Result<()> {
        let file: File = File::create(path)?;
        let writer: BufWriter<File> = BufWriter::new(file);
        
        bincode::serialize_into(writer, self)
            .map_err(|e: Box<bincode::ErrorKind>| io::Error::new(io::ErrorKind::Other, format!("Failed to serialize graph: {}", e)))
    }

    /// Deserialize the graph from a binary file into memory
    /// Example: Load and analyze
    /// let loaded_graph = PanGenomeGraph::load_from_file("virus_pangenome.bin")?;
    /// println!("Loaded graph with {} nodes!", loaded_graph.graph.node_count());
    /// 
    pub fn load_from_file(path: &str) -> io::Result<Self> {
        let file: File = File::open(path)?;
        let reader: BufReader<File> = BufReader::new(file);
        
        bincode::deserialize_from(reader)
            .map_err(|e: Box<bincode::ErrorKind>| io::Error::new(io::ErrorKind::Other, format!("Failed to deserialize graph: {}", e)))
    }

    /// Usage:
    /// fn main() {
    /// let k = 31; // Standard k-mer size
    ///
    /// match PanGenomeGraph::from_fastas("reference.fasta", "assemblies.fasta", k) {
    ///     Ok(graph) => {
    ///         println!(
    ///             "Successfully built Pan-Genome Graph!\nNodes: {}\nEdges: {}",
    ///             graph.graph.node_count(),
    ///             graph.graph.edge_count()
    ///         );
    ///     }
    ///     Err(e) => {
    ///         eprintln!("Failed to build graph due to an IO error: {}", e);
    ///     }
    /// }
    /// 
    pub fn from_fastas(reference_path: &str, assemblies_path: &str, k: usize) -> io::Result<Self> {
        let mut graph: PanGenomeGraph = Self::new(k);

        // Ingest the reference backbone
        let ref_reader: FastaReader<crate::readers::ReaderType> = FastaReader::from_file(reference_path)?;
        // Assume the first sequence in the reference file is the backbone
        let ref_record: crate::fasta::FastaRecord = ref_reader.into_iter().next().expect("Reference file is empty")?;
        let ref_bytes: &[u8] = ref_record.seq.as_bytes();
        // Add the reference to the graph
        graph.add_sequence(std::str::from_utf8(ref_bytes).unwrap());
        
        // Create the orientor using a small k-mer for flexible mapping (e.g., 15)
        let orientor: StrandOrientor = StrandOrientor::new(ref_bytes, 15);

        // Add assemblies to the graph
        let assembly_reader: FastaReader<crate::readers::ReaderType> = FastaReader::from_file(assemblies_path)?;
        for result in assembly_reader {
            let record: crate::fasta::FastaRecord = result?;
            if !record.is_empty() {
                // Instantly orient the sequence to match the reference strand
                let oriented_bytes: Vec<u8> = orientor.orient(record.seq.as_bytes());
                
                // Convert back to string and add to graph
                let oriented_str: &str = std::str::from_utf8(&oriented_bytes)
                    .expect("Invalid UTF-8 after orientation");
                    
                graph.add_sequence(oriented_str);
            }
        }

        Ok(graph)
    }
}

pub struct StrandOrientor {
    /// A set of strictly directional (forward) 2-bit encoded k-mers from the reference
    reference_kmers: HashSet<u64>,
    k: usize,
}

impl StrandOrientor {
    /// Initializes the orientor using the reference backbone
    pub fn new(reference: &[u8], k: usize) -> Self {
        let mut reference_kmers = HashSet::new();
        
        for kmer in reference.to_kmers(k) {
            // Use encode_to_u64(), NOT canonical_u64(), for orientation
            reference_kmers.insert(kmer.encode_to_u64());
        }
        
        Self { reference_kmers, k }
    }

    /// Test a sequence and return it correctly oriented to the reference strand
    pub fn orient(&self, sequence: &[u8]) -> Vec<u8> {
        let mut fwd_hits: i32 = 0;
        let mut rc_hits: i32  = 0;

        let rc_seq: Vec<u8> = sequence.reverse_complement();

        // 1. Count hits for the forward sequence
        for kmer in sequence.to_kmers(self.k) {
            if self.reference_kmers.contains(&kmer.encode_to_u64()) {
                fwd_hits += 1;
            }
        }

        // 2. Count hits for the reverse complement sequence
        for kmer in rc_seq.to_kmers(self.k) {
            if self.reference_kmers.contains(&kmer.encode_to_u64()) {
                rc_hits += 1;
            }
        }

        // 3. Return whichever sequence maps better to the reference strand
        if rc_hits > fwd_hits {
            rc_seq
        } else {
            sequence.to_vec()
        }
    }
}