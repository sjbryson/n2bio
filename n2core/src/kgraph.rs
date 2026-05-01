//! n2core/src/kgraph.rs

use std::fs::File;
use std::io::{self, BufReader, BufWriter};
use serde::{Deserialize, Serialize};
use petgraph::graph::{DiGraph, NodeIndex};
use std::collections::HashMap;
use crate::sequence::{DnaSequence, Kmer};
use crate::fasta::FastaReader;


/// Data stored in each node
#[derive(Serialize, Deserialize)]
pub struct KmerNode {
    pub canonical_seq: Vec<u8>,
    pub frequency: u32,
    // Use HashSet<String> to track genomes associated with each kmer?
}
#[derive(Serialize, Deserialize)]
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

    /// Serialize the graph to a binary file
    /// Example: Build and save
    /// let graph = PanGenomeGraph::from_fastas("ref.fa", "assemblies.fa", 31)?;
    /// graph.save_to_file("virus_pangenome.bin")?;
    /// 
    pub fn save_to_file(&self, path: &str) -> io::Result<()> {
        let file = File::create(path)?;
        let writer  = BufWriter::new(file);
        
        bincode::serialize_into(writer, self)
            .map_err(|e| io::Error::new(io::ErrorKind::Other, format!("Failed to serialize graph: {}", e)))
    }

    /// Deserialize the graph from a binary file into memory
    /// Example: Load and analyze
    /// let loaded_graph = PanGenomeGraph::load_from_file("virus_pangenome.bin")?;
    /// println!("Loaded graph with {} nodes!", loaded_graph.graph.node_count());
    /// 
    pub fn load_from_file(path: &str) -> io::Result<Self> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        
        bincode::deserialize_from(reader)
            .map_err(|e| io::Error::new(io::ErrorKind::Other, format!("Failed to deserialize graph: {}", e)))
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

