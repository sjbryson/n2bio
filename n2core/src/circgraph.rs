//! n2core/src/circgraph.rs

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::collections::HashMap;
use serde::{Deserialize, Serialize};
use petgraph::graph::{DiGraph, NodeIndex};
use petgraph::visit::EdgeRef;
use crate::sequence::DnaSequence;
use crate::kmer::{KmerEncoding, StrandOrientor};
use crate::fasta::FastaReader;

// ============================================================================
// NEEDLEMAN-WUNSCH ALIGNMENT ENGINE
// ============================================================================

#[derive(Debug, Clone, PartialEq)]
pub enum EditOp {
    /// Both sequences have the exact same base
    Match,
    /// The assembly has a different base (Substitution)
    Substitution(u8),
    /// The assembly has an extra base not in the reference (Insertion)
    Insertion(u8),
    /// The assembly is missing a base present in the reference (Deletion)
    Deletion,
}

/// A cache-friendly, 1D flattened Needleman-Wunsch implementation 
/// specifically designed for small gap filling between anchors.
pub fn align_gaps(reference: &[u8], assembly: &[u8]) -> Vec<EditOp> {
    let r_len = reference.len();
    let a_len = assembly.len();

    let match_score: i32 = 1;
    let mismatch_score: i32 = -1;
    let gap_penalty: i32 = -1;

    let cols = a_len + 1;
    let rows = r_len + 1;

    let mut dp = vec![0i32; rows * cols];
    let mut tb = vec![0u8; rows * cols]; // 0=Diag, 1=Up, 2=Left

    for i in 1..rows {
        dp[i * cols] = i as i32 * gap_penalty;
        tb[i * cols] = 1; 
    }
    for j in 1..cols {
        dp[j] = j as i32 * gap_penalty;
        tb[j] = 2; 
    }

    for i in 1..rows {
        for j in 1..cols {
            let diag_score = dp[(i - 1) * cols + (j - 1)] + 
                if reference[i - 1] == assembly[j - 1] { match_score } else { mismatch_score };
            let up_score = dp[(i - 1) * cols + j] + gap_penalty;
            let left_score = dp[i * cols + (j - 1)] + gap_penalty;

            let max_score = diag_score.max(up_score).max(left_score);
            dp[i * cols + j] = max_score;

            if max_score == diag_score {
                tb[i * cols + j] = 0;
            } else if max_score == up_score {
                tb[i * cols + j] = 1;
            } else {
                tb[i * cols + j] = 2;
            }
        }
    }

    let mut i = r_len;
    let mut j = a_len;
    let mut ops = Vec::with_capacity(r_len.max(a_len));

    while i > 0 || j > 0 {
        match tb[i * cols + j] {
            0 => {
                if reference[i - 1] == assembly[j - 1] {
                    ops.push(EditOp::Match);
                } else {
                    ops.push(EditOp::Substitution(assembly[j - 1]));
                }
                i -= 1;
                j -= 1;
            }
            1 => {
                ops.push(EditOp::Deletion);
                i -= 1;
            }
            2 => {
                ops.push(EditOp::Insertion(assembly[j - 1]));
                j -= 1;
            }
            _ => unreachable!(),
        }
    }

    ops.reverse();
    ops
}

// ============================================================================
// CIRCULAR VARIATION GRAPH ARCHITECTURE
// ============================================================================

#[derive(Serialize, Deserialize, Debug, Clone)]
pub enum NodeType {
    Nucleotide {
        base: u8,
        ref_pos: Option<usize>,
    },
    Feature {
        name: String,
        kind: String,
    },
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct PangenomeNode {
    pub data: NodeType,
    pub frequency: u32,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub enum EdgeType {
    Sequence,
    AnnotationStart,
    AnnotationStop,
    VariantBranch,
}

#[derive(Serialize, Deserialize)]
pub struct CircularVariationGraph {
    pub graph: DiGraph<PangenomeNode, EdgeType>,
    pub reference_backbone: Vec<NodeIndex>,
    pub anchors: HashMap<u64, NodeIndex>,
    pub k: usize,
}

impl CircularVariationGraph {
    pub fn new(k: usize) -> Self {
        assert!(k <= 32, "ERROR: K-mer size must be <= 32.");
        assert!(k > 0,   "ERROR: K-mer size must be > 0.");

        Self {
            graph: DiGraph::new(),
            reference_backbone: Vec::new(),
            anchors: HashMap::new(),
            k,
        }
    }

     /// Try to create a new graph, return Error if k is out of bounds
    pub fn try_new(k: usize) -> Result<Self, String> {
        if k > 32 {
            return Err(format!("ERROR: K-mer size (k={}) must be <= 32.", k));
        }
        if k == 0 {
            return Err("ERROR: K-mer size must be > 0.".to_string());
        }

        Ok(Self {
            graph: DiGraph::new(),
            reference_backbone: Vec::new(),
            anchors: HashMap::new(),
            k,
        })
    }

    /// Initializes the circular graph with the reference backbone
    pub fn add_reference_circular(&mut self, seq: &[u8]) {
        let l: usize = seq.len();
        self.reference_backbone.reserve(l);
        let mut prev_node: Option<NodeIndex> = None;

        // 1. Build the Linear Backbone
        for (pos, &base) in seq.iter().enumerate() {
            let node = PangenomeNode {
                data: NodeType::Nucleotide { base, ref_pos: Some(pos) },
                frequency: 1,
            };
            
            let current_node: NodeIndex = self.graph.add_node(node);
            self.reference_backbone.push(current_node);

            if let Some(prev) = prev_node {
                self.graph.add_edge(prev, current_node, EdgeType::Sequence);
            }
            prev_node = Some(current_node);
        }

        // 2. CLOSE THE CIRCLE!
        if l > 0 {
            let first_node: NodeIndex = self.reference_backbone[0];
            let last_node: NodeIndex = *self.reference_backbone.last().unwrap();
            self.graph.add_edge(last_node, first_node, EdgeType::Sequence);
        }

        // 3. Find Unique Anchors (Including origin-spanning k-mers)
        let mut extended_seq: Vec<u8> = seq.to_vec();
        extended_seq.extend_from_slice(&seq[..self.k - 1]);

        let mut kmer_counts: HashMap<u64, (usize, i32)> = HashMap::new();
        for (i, kmer) in extended_seq.to_kmers(self.k).enumerate() {
            let encoded: u64 = kmer.encode_to_u64();
            let pos: usize = i % l; 
            let entry: &mut (usize, i32) = kmer_counts.entry(encoded).or_insert((pos, 0));
            entry.1 += 1;
        }

        for (encoded, (pos, count)) in kmer_counts {
            if count == 1 {
                self.anchors.insert(encoded, self.reference_backbone[pos]);
            }
        }
    }

    /// Extracts reference bases by walking the circular graph
    fn get_ref_bytes_cyclic(&self, start: usize, length: usize) -> Vec<u8> {
        let l: usize = self.reference_backbone.len();
        let mut bytes: Vec<u8> = Vec::with_capacity(length);
        
        for i in 0..length {
            let r_idx: usize = (start + i) % l;
            if let NodeType::Nucleotide { base, .. } = self.graph[self.reference_backbone[r_idx]].data {
                bytes.push(base);
            }
        }
        bytes
    }

    /// Finds an existing variant branch or creates a new one
    fn find_or_create_variant(&mut self, source: Option<NodeIndex>, base: u8, ref_pos: Option<usize>) -> NodeIndex {
        if let Some(src) = source {
            for edge in self.graph.edges(src) {
                let target: NodeIndex = edge.target();
                if let NodeType::Nucleotide { base: t_base, ref_pos: t_ref_pos } = self.graph[target].data {
                    if t_base == base && t_ref_pos == ref_pos {
                        self.graph[target].frequency += 1;
                        return target;
                    }
                }
            }
            
            let new_node: PangenomeNode = PangenomeNode { data: NodeType::Nucleotide { base, ref_pos }, frequency: 1 };
            let new_idx: NodeIndex = self.graph.add_node(new_node);
            self.graph.add_edge(src, new_idx, EdgeType::VariantBranch);
            new_idx
        } else {
            let new_node: PangenomeNode = PangenomeNode { data: NodeType::Nucleotide { base, ref_pos }, frequency: 1 };
            self.graph.add_node(new_node)
        }
    }

    /// Applies EditOps to the circular graph using modulo traversal
    fn apply_alignment_cyclic(&mut self, ops: Vec<EditOp>, mut curr_r_unrolled: usize, mut current_tail: Option<NodeIndex>) -> Option<NodeIndex> {
        let l: usize = self.reference_backbone.len();
        for op in ops {
            match op {
                EditOp::Match => {
                    let next_node: NodeIndex = self.reference_backbone[curr_r_unrolled % l];
                    self.graph[next_node].frequency += 1;
                    if let Some(tail) = current_tail {
                        if self.graph.find_edge(tail, next_node).is_none() {
                            self.graph.add_edge(tail, next_node, EdgeType::VariantBranch);
                        }
                    }
                    current_tail = Some(next_node);
                    curr_r_unrolled += 1;
                }
                EditOp::Substitution(b) => {
                    let new_tail: NodeIndex = self.find_or_create_variant(current_tail, b, Some(curr_r_unrolled % l));
                    current_tail = Some(new_tail);
                    curr_r_unrolled += 1;
                }
                EditOp::Insertion(b) => {
                    let new_tail: NodeIndex = self.find_or_create_variant(current_tail, b, None);
                    current_tail = Some(new_tail);
                }
                EditOp::Deletion => {
                    curr_r_unrolled += 1;
                }
            }
        }
        current_tail
    }

    /// Aligns and merges a circular sequence with origin-spanning and unrolled chaining
    pub fn add_circular_sequence(&mut self, sequence: &[u8]) {
        let l: usize = self.reference_backbone.len();
        if sequence.len() < self.k || l == 0 { return; }

        // --- Phase 1: Seeding ---
        let mut found_anchors: Vec<(usize, usize)> = Vec::new(); 
        for (asm_pos, kmer) in sequence.to_kmers(self.k).enumerate() {
            if let Some(&node_idx) = self.anchors.get(&kmer.encode_to_u64()) {
                if let NodeType::Nucleotide { ref_pos: Some(r_pos), .. } = self.graph[node_idx].data {
                    found_anchors.push((asm_pos, r_pos));
                }
            }
        }
        
        if found_anchors.is_empty() { return; }
        found_anchors.sort_by_key(|&(a, _)| a);

        // --- Phase 2: Circular Unrolling & Chaining ---
        let mut chain: Vec<(usize, usize)> = Vec::new();
        let (mut last_a, mut last_unrolled_r) = found_anchors[0];
        chain.push((last_a, last_unrolled_r));

        for &(a, r) in found_anchors.iter().skip(1) {
            let a_dist: usize = a - last_a;
            let r_dist: usize = (r + l - (last_unrolled_r % l)) % l;

            if a_dist >= self.k && r_dist >= self.k {
                let unrolled_r: usize = last_unrolled_r + r_dist;
                chain.push((a, unrolled_r));
                last_a = a;
                last_unrolled_r = unrolled_r;
            }
        }

        let mut current_tail_node: Option<NodeIndex> = None;

        // --- Phase 3.1: The 5' Tail (Wrapping Backwards) ---
        let (first_a, first_unrolled_r) = chain[0];
        if first_a > 0 {
            let asm_tail: &[u8] = &sequence[0..first_a];
            let ref_start_unrolled: usize = first_unrolled_r.saturating_sub(first_a);
            let ref_gap_len: usize = first_unrolled_r - ref_start_unrolled;
            let ref_start_modulo: usize = ref_start_unrolled % l;
            
            let ref_tail: Vec<u8> = self.get_ref_bytes_cyclic(ref_start_modulo, ref_gap_len);
            let ops: Vec<EditOp> = align_gaps(&ref_tail, asm_tail);
            
            let starting_node: Option<NodeIndex> = Some(self.reference_backbone[(ref_start_modulo + l - 1) % l]);
            current_tail_node = self.apply_alignment_cyclic(ops, ref_start_unrolled, starting_node);
        }

        // --- Phase 3.2: The Internal Anchors and Gaps ---
        for i in 0..chain.len() {
            let (a, unrolled_r) = chain[i];
            let true_r: usize = unrolled_r % l;

            if i > 0 {
                let (prev_a, prev_unrolled_r) = chain[i - 1];
                let asm_gap: &[u8] = &sequence[prev_a + self.k .. a];
                
                let ref_start_unrolled: usize = prev_unrolled_r + self.k;
                let ref_gap_len: usize = unrolled_r.saturating_sub(ref_start_unrolled);
                
                let ref_gap: Vec<u8> = self.get_ref_bytes_cyclic(ref_start_unrolled % l, ref_gap_len);
                let ops: Vec<EditOp> = align_gaps(&ref_gap, asm_gap);
                
                current_tail_node = self.apply_alignment_cyclic(ops, ref_start_unrolled, current_tail_node);
            }

            let anchor_start_node = self.reference_backbone[true_r];
            if let Some(tail) = current_tail_node {
                if self.graph.find_edge(tail, anchor_start_node).is_none() {
                    self.graph.add_edge(tail, anchor_start_node, EdgeType::VariantBranch);
                }
            }

            for k_offset in 0..self.k {
                let node_idx: NodeIndex = self.reference_backbone[(true_r + k_offset) % l];
                self.graph[node_idx].frequency += 1;
                current_tail_node = Some(node_idx);
            }
        }

        // --- Phase 3.3: The 3' Tail (Wrapping Forwards) ---
        let (last_a, last_unrolled_r) = chain.last().unwrap();
        let asm_tail_start: usize = last_a + self.k;
        
        if asm_tail_start < sequence.len() {
            let asm_tail: &[u8] = &sequence[asm_tail_start..];
            let ref_start_unrolled: usize = last_unrolled_r + self.k;
            
            // For a circular graph, we can always extract enough reference bases
            let ref_tail: Vec<u8> = self.get_ref_bytes_cyclic(ref_start_unrolled % l, asm_tail.len());
            let ops: Vec<EditOp> = align_gaps(&ref_tail, asm_tail);
            
            self.apply_alignment_cyclic(ops, ref_start_unrolled, current_tail_node);
        }
    }

    pub fn from_fastas(reference_path: &str, assemblies_path: &str, k: usize) -> io::Result<Self> {
        // Initialize the new circular graph safely with try_new
        let mut graph: CircularVariationGraph = Self::try_new(k).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        // Ingest the reference backbone
        let ref_reader: FastaReader<crate::readers::ReaderType> = FastaReader::from_file(reference_path)?;
        
        // Assume the first sequence in the reference file is the backbone
        let ref_record: crate::fasta::FastaRecord = ref_reader.into_iter().next().expect("Reference file is empty")?;
        let ref_bytes: &[u8] = ref_record.seq.as_bytes();
        
        // Add the reference to the graph using the circular initialization
        graph.add_reference_circular(ref_bytes);
        
        // Create the orientor using a small k-mer for flexible mapping (e.g., 15)
        let orientor: StrandOrientor = StrandOrientor::new(ref_bytes, 15);

        // Add assemblies to the graph
        let assembly_reader: FastaReader<crate::readers::ReaderType> = FastaReader::from_file(assemblies_path)?;
        
        for result in assembly_reader {
            let record: crate::fasta::FastaRecord = result?;
            if !record.is_empty() {
                
                // Orient the sequence to match the reference strand
                let oriented_bytes: Vec<u8> = orientor.orient(record.seq.as_bytes());
                
                // Add the aligned assembly directly as bytes using the cyclic method
                graph.add_circular_sequence(&oriented_bytes);
            }
        }

        Ok(graph)
    }


    /// Exports to Cytoscape-friendly GraphML
    pub fn export_to_graphml(&self, path: &str) -> io::Result<()> {
        let file: File = File::create(path)?;
        let mut writer: BufWriter<File> = BufWriter::new(file);

        writeln!(writer, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>")?;
        writeln!(writer, "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\">")?;
        
        writeln!(writer, "  <key id=\"v_type\" for=\"node\" attr.name=\"type\" attr.type=\"string\"/>")?;
        writeln!(writer, "  <key id=\"v_label\" for=\"node\" attr.name=\"label\" attr.type=\"string\"/>")?;
        writeln!(writer, "  <key id=\"v_freq\" for=\"node\" attr.name=\"frequency\" attr.type=\"int\"/>")?;
        writeln!(writer, "  <key id=\"v_pos\" for=\"node\" attr.name=\"ref_pos\" attr.type=\"int\"/>")?;
        writeln!(writer, "  <key id=\"e_type\" for=\"edge\" attr.name=\"type\" attr.type=\"string\"/>")?;

        writeln!(writer, "  <graph id=\"circular_variation_graph\" edgedefault=\"directed\">")?;

        for node_idx in self.graph.node_indices() {
            let node: &PangenomeNode = &self.graph[node_idx];
            writeln!(writer, "    <node id=\"n{}\">", node_idx.index())?;
            writeln!(writer, "      <data key=\"v_freq\">{}</data>", node.frequency)?;
            
            match &node.data {
                NodeType::Nucleotide { base, ref_pos } => {
                    writeln!(writer, "      <data key=\"v_type\">Nucleotide</data>")?;
                    writeln!(writer, "      <data key=\"v_label\">{}</data>", *base as char)?;
                    if let Some(pos) = ref_pos {
                        writeln!(writer, "      <data key=\"v_pos\">{}</data>", pos)?;
                    }
                }
                NodeType::Feature { name, kind } => {
                    writeln!(writer, "      <data key=\"v_type\">Feature_{}</data>", kind)?;
                    writeln!(writer, "      <data key=\"v_label\">{}</data>", name)?;
                }
            }
            writeln!(writer, "    </node>")?;
        }

        for edge in self.graph.edge_references() {
            let source: usize = edge.source().index();
            let target: usize = edge.target().index();
            let e_type: &str = match edge.weight() {
                EdgeType::Sequence        => "Sequence",
                EdgeType::AnnotationStart => "Start",
                EdgeType::AnnotationStop  => "Stop",
                EdgeType::VariantBranch   => "Variant",
            };

            writeln!(writer, "    <edge source=\"n{}\" target=\"n{}\">", source, target)?;
            writeln!(writer, "      <data key=\"e_type\">{}</data>", e_type)?;
            writeln!(writer, "    </edge>")?;
        }

        writeln!(writer, "  </graph>\n</graphml>")?;
        Ok(())
    }
}
