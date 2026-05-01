//! n2core/src/kgraph.rs

use std::io;
use petgraph::graph::{DiGraph, NodeIndex};
use std::collections::HashMap;
use crate::sequence::{DnaSequence, Kmer};
use crate::fasta::FastaReader;


/// Data stored in each node
pub struct KmerNode {
    pub canonical_seq: Vec<u8>,
    pub frequency: u32,
    // Use HashSet<String> to track genomes associated with each kmer?
}

pub struct PanGenomeGraph {
    /// Petgraph storing nodes and edges
    pub graph: DiGraph<KmerNode, u32>,
    /// Lookup index mapping canonical k-mers to nodes
    pub node_map: HashMap<Vec<u8>, NodeIndex>,
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
            let canonical: Vec<u8> = kmer.canonical();

            // 1. Get the existing node, or create a new one
            let current_node: NodeIndex = if let Some(&idx) = self.node_map.get(&canonical) {
                self.graph[idx].frequency += 1;
                idx
            } else {
                let new_node: KmerNode = KmerNode {
                    canonical_seq: canonical.clone(),
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

    ///Usage:
    ///fn main() {
    ///let k = 31; // Standard k-mer size
    ///
    ///match PanGenomeGraph::from_fastas("reference.fasta", "assemblies.fasta", k) {
    ///    Ok(graph) => {
    ///        println!(
    ///            "Successfully built Pan-Genome Graph!\nNodes: {}\nEdges: {}",
    ///            graph.graph.node_count(),
    ///            graph.graph.edge_count()
    ///        );
    ///    }
    ///    Err(e) => {
    ///        eprintln!("Failed to build graph due to an IO error: {}", e);
    ///    }
    ///}
    pub fn from_fastas(reference_path: &str, assemblies_path: &str, k: usize) -> io::Result<Self> {
        let mut graph: PanGenomeGraph = Self::new(k);

        // 1. Ingest the reference backbone
        let ref_reader: FastaReader<crate::readers::ReaderType> = FastaReader::from_file(reference_path)?;
        for result in ref_reader {
            let record: crate::fasta::FastaRecord = result?;
            if !record.is_empty() {
                graph.add_sequence(&record.seq);
            }
        }

        // 2. Layer the Assemblies on top
        let assembly_reader: FastaReader<crate::readers::ReaderType> = FastaReader::from_file(assemblies_path)?;
        for result in assembly_reader {
            let record: crate::fasta::FastaRecord = result?;
            if !record.is_empty() {
                graph.add_sequence(&record.seq);
            }
        }

        Ok(graph)
    }
}

