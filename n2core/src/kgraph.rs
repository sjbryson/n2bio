//! n2core/src/kgraph.rs

use std::fs::File;
use std::io::{self, BufReader, BufWriter, Write};
use serde::{Deserialize, Serialize};
use petgraph::graph::{DiGraph, NodeIndex};
use petgraph::visit::EdgeRef;
use std::collections::HashMap;
use crate::sequence::DnaSequence;
use crate::kmer::{KmerEncoding, StrandOrientor};
use crate::fasta::FastaReader;


/// Data stored in each node
#[derive(Serialize, Deserialize)]
pub struct KmerNode {
    pub kmer_u64: u64,
    pub frequency: u32,
    // Use HashSet<String> to track genomes associated with each kmer?
}

#[derive(Serialize, Deserialize)]
pub struct CanonicalKmerNode {
    pub canonical_u64: u64,
    pub frequency: u32,
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
        // Enforce k-mer size for 2-bit encoding
        assert!(k <= 32, "ERROR: K-mer size (k={}) exceeds the maximum of 32 for u64 encoding.", k);
        assert!(k > 0,   "ERROR: K-mer size must be greater than 0.");

        Self {
            graph: DiGraph::new(),
            node_map: HashMap::new(),
            k,
        }
    }

    /// Try to create a new graph, return Error if k is out of bounds
    pub fn try_new(k: usize) -> Result<Self, String> {
        if k > 32 {
            return Err(format!("K-mer size (k={}) exceeds the maximum of 32 for u64 encoding.", k));
        }
        if k == 0 {
            return Err("K-mer size must be greater than 0.".to_string());
        }

        Ok(Self {
            graph: DiGraph::new(),
            node_map: HashMap::new(),
            k,
        })
    }

    /// Incorporates a new sequence in the graph
    pub fn add_sequence(&mut self, seq: &str) {
        let mut prev_node: Option<NodeIndex> = None;

        // Calling .as_bytes(), uses the &[u8] implementation of DnaSequence
        for kmer in seq.as_bytes().to_kmers(self.k) {
            let kmer_encoded: u64 = kmer.encode_to_u64();

            // 1. Get the existing node, or create a new one
            let current_node: NodeIndex = if let Some(&idx) = self.node_map.get(&kmer_encoded) {
                self.graph[idx].frequency += 1;
                idx
            } else {
                let new_node: KmerNode = KmerNode {
                    kmer_u64: kmer_encoded.clone(),
                    frequency: 1,
                };
                let idx: NodeIndex = self.graph.add_node(new_node);
                self.node_map.insert(kmer_encoded, idx);
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

    /// Exports the graph to a Cytoscape-compatible GraphML file
    pub fn export_to_graphml(&self, path: &str) -> io::Result<()> {
        let file: File = File::create(path)?;
        let mut writer: BufWriter<File> = BufWriter::new(file);

        // 1. Write the standard GraphML XML headers
        writeln!(writer, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>")?;
        writeln!(writer, "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\"")?;
        writeln!(writer, "         xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"")?;
        writeln!(writer, "         xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">")?;
        
        // 2. Define Node Data Keys (Sequence and Frequency)
        writeln!(writer, "  <key id=\"v_seq\" for=\"node\" attr.name=\"sequence\" attr.type=\"string\"/>")?;
        writeln!(writer, "  <key id=\"v_freq\" for=\"node\" attr.name=\"frequency\" attr.type=\"int\"/>")?;
        
        // 3. Define Edge Data Keys (Weight / Transition Frequency)
        writeln!(writer, "  <key id=\"e_weight\" for=\"edge\" attr.name=\"weight\" attr.type=\"int\"/>")?;

        // 4. Open the Directed Graph Block
        writeln!(writer, "  <graph id=\"pangenome\" edgedefault=\"directed\">")?;

        // 5. Iterate and write Nodes
        for node_idx in self.graph.node_indices() {
            let node_data = &self.graph[node_idx];
            
            // Decode the u64 back to the kmer sequence
            let seq_bytes: Vec<u8> = <[u8]>::decode_from_u64(node_data.kmer_u64, self.k);
            let seq_str: &str = std::str::from_utf8(&seq_bytes)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

            // Write the XML block for the node using petgraph's internal index as the ID
            writeln!(writer, "    <node id=\"n{}\">", node_idx.index())?;
            writeln!(writer, "      <data key=\"v_seq\">{}</data>", seq_str)?;
            writeln!(writer, "      <data key=\"v_freq\">{}</data>", node_data.frequency)?;
            writeln!(writer, "    </node>")?;
        }

        // 6. Iterate and write Edges
        for edge in self.graph.edge_references() {
            let source_idx: usize = edge.source().index();
            let target_idx: usize = edge.target().index();
            let weight: &u32      = edge.weight(); // This is the u32 frequency payload

            // Edge IDs are optional in GraphML, so we just declare source and target
            writeln!(writer, "    <edge source=\"n{}\" target=\"n{}\">", source_idx, target_idx)?;
            writeln!(writer, "      <data key=\"e_weight\">{}</data>", weight)?;
            writeln!(writer, "    </edge>")?;
        }

        // 7. Close the GraphML tags
        writeln!(writer, "  </graph>")?;
        writeln!(writer, "</graphml>")?;

        Ok(())
    }

    /// Usage:
    ///
    ///     let k = 31; // Standard k-mer size
    ///
    ///     match PanGenomeGraph::from_fastas("reference.fasta", "assemblies.fasta", k) {
    ///         Ok(graph) => {
    ///             println!(
    ///                 "Successfully built Pan-Genome Graph!\nNodes: {}\nEdges: {}",
    ///                  graph.graph.node_count(),
    ///                  graph.graph.edge_count()
    ///             );
    ///         }
    ///         Err(e) => {
    ///             eprintln!("Failed to build graph due to an IO error: {}", e);
    ///         }
    /// 
    pub fn from_fastas(reference_path: &str, assemblies_path: &str, k: usize) -> io::Result<Self> {
        let mut graph: PanGenomeGraph = Self::try_new(k).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

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
                
                // Orient the sequence to match the reference strand
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
