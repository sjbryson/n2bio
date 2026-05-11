//! n2core/src/alignmentgraph.rs
//!

use std::collections::{HashMap};
use std::fs::File;
use std::io::{self, BufWriter, Write};
use thiserror::Error;

use crate::kmer::{ KmerEncoding, StrandOrientor, KmerError };
use crate::sequence::DnaSequence;
use crate::align::{ EditOp, AlignmentScores, align_seqs, align_poa };

// ============================================================================
// Errors
// ============================================================================

#[derive(Error, Debug)]
pub enum AlignmentGraphError {
    #[error("error")]
    AlignmentGraphError(),
}

// ============================================================================
// Edge Definition
// ============================================================================

/// A bidirectional edge connecting two nodes in the sequence graph.
#[derive(Debug, Clone)]
pub struct Edge {
    pub from_node:   usize,
    pub from_strand: Strand,
    pub to_node:     usize,
    pub to_strand:   Strand,
    pub edge_type:   EdgeType,
}

/// Categorizes the biological meaning of an edge.
#[derive(Debug, Clone)]
pub enum EdgeType {
    Backbone,  // The original reference genome sequence
    Match,     // Query genome matches the reference
    Mismatch,  // Query genome has a SNP/substitution
    Insertion, // Query genome inserts sequence not in reference
    Deletion,  // Query genome skips over reference sequence
}

// ============================================================================
// Alignment graph
// ============================================================================

#[derive(Debug, Clone)]
pub struct AlignmentGraph {
    pub positions:   HashMap<usize, ReferencePosition>,
    pub anchors:     HashMap<u64, ReferenceAnchor>,
    pub edges_out:   HashMap<usize, Vec<Edge>>,
    pub edges_in:    HashMap<usize, Vec<Edge>>,
    pub is_circular: bool,
    //pub anchor_len -> kmer size 
    next_node_id:    usize,
}

impl AlignmentGraph {
    pub fn new(circular: bool) -> Result<Self, AlignmentGraphError> {
        Ok(AlignmentGraph {
            positions:    HashMap::new(),
            anchors:      HashMap::new(),
            edges_out:    HashMap::new(),
            edges_in:     HashMap::new(),
            is_circular:  circular,
            next_node_id: 0,
        })
    }
}

// ============================================================================
// Alignment graph - Structure & Orientation
// ============================================================================

impl AlignmentGraph {
    pub fn add_reference(&mut self, sequence: &[u8]) {
        if !sequence.is_empty() {
            // 1. Create all ReferencePosition nodes
            for (pos, &nuc) in sequence.iter().enumerate() {
                let node_id: usize = self.add_position_node(pos, nuc);
                
                // 2. Create directed backbone edges
                if pos > 0 {
                    let prev_node_id = node_id - 1;
                    self.add_edge(
                        prev_node_id, Strand::Forward, 
                        node_id, Strand::Forward, 
                        EdgeType::Backbone
                    );
                }
            }
            // 3. Handle optional circularization
            if self.is_circular && sequence.len() > 1 {
                let first_node_id: usize = 0;
                let last_node_id:  usize = sequence.len() - 1;
                self.add_edge(
                    last_node_id, Strand::Forward, 
                    first_node_id, Strand::Forward, 
                    EdgeType::Backbone
                );
            }
        }
    }

    pub fn add_genome(
        &mut self,
        query_seq: &[u8],
        k: usize,
        orientor: &StrandOrientor,
        scores: &AlignmentScores,
    ) -> Result<(), Box<dyn std::error::Error>> {
        
        let (oriented_query, matches) = self.map_query_anchors(query_seq, orientor, k)?;
        let segments = self.get_unaligned_segments(&matches, k, oriented_query.len());
        
        for seg in segments {
            match seg.segment_type {
                SegmentType::Bounded { ref_node_start, ref_node_end } => {
                    let query_slice = &oriented_query[seg.query_start..seg.query_end];
                    
                    // 1. Extract the DAG bubble zone between the anchors
                    let sorted_nodes = self.get_topological_subgraph(ref_node_start, ref_node_end);
                    let (node_bases, predecessors) = self.build_poa_inputs(&sorted_nodes);

                    // 2. Align against the DAG using POA
                    let (ops, poa_path) = align_poa(query_slice, &node_bases, &predecessors, scores);
                    
                    // 3. Weave it back in
                    let anchor_start_node = ref_node_start.saturating_sub(1);
                    self.weave_poa_alignment(&ops, &poa_path, &sorted_nodes, Some(anchor_start_node), Some(ref_node_end));
                }
                
                // For extending flanks on linear genomes, we can safely fall back to standard 
                // Needleman-Wunsch, as structural variation generally lives between anchors.
                SegmentType::ExtendBackward { ref_node_end } => {
                    let query_slice = &oriented_query[seg.query_start..seg.query_end];
                    let ref_path: Vec<usize> = (0..ref_node_end).collect();
                    let ref_slice: Vec<u8> = ref_path.iter().filter_map(|&id| self.positions.get(&id).map(|n| n.nucleotide)).collect();

                    let rev_query: Vec<u8> = query_slice.iter().rev().copied().collect();
                    let rev_ref: Vec<u8> = ref_slice.iter().rev().copied().collect();

                    let mut ops = align_seqs(&rev_ref, &rev_query, scores);
                    ops.reverse(); 

                    self.weave_alignment(&ops, &ref_path, None, Some(ref_node_end));
                }
                
                SegmentType::ExtendForward { ref_node_start } => {
                    let query_slice = &oriented_query[seg.query_start..seg.query_end];
                    // Approximate extending to the end of the original backbone
                    let last_node = self.positions.values().map(|n| n.id).max().unwrap_or(0);
                    let ref_path: Vec<usize> = (ref_node_start..=last_node).collect();
                    let ref_slice: Vec<u8> = ref_path.iter().filter_map(|&id| self.positions.get(&id).map(|n| n.nucleotide)).collect();

                    let ops = align_seqs(&ref_slice, query_slice, scores);
                    let anchor_start_node = ref_node_start.saturating_sub(1);

                    self.weave_alignment(&ops, &ref_path, Some(anchor_start_node), None);
                }
            }
        }
        Ok(())
    }

    // Helper to add nodes (updated to initialize edge lists)
    fn add_position_node(&mut self, position: usize, nucleotide: u8) -> usize {
        let id: usize = self.next_node_id;
        self.positions.insert(id, ReferencePosition { id, position, nucleotide });
        self.edges_out.insert(id, Vec::new());
        self.edges_in.insert(id, Vec::new());
        self.next_node_id += 1;
        id
    }

    /// Adds a bidirected edge. It registers the outgoing edge for the 'from' node
    /// and the incoming edge for the 'to' node.
    pub fn add_edge(&mut self, from_node: usize, from_strand: Strand, to_node: usize, to_strand: Strand, edge_type: EdgeType) {
        let edge: Edge = Edge { from_node, from_strand, to_node, to_strand, edge_type };
        
        // Outgoing from the perspective of the 'from' node
        if let Some(out_list) = self.edges_out.get_mut(&from_node) {
            out_list.push(edge.clone());
        }
        
        // Incoming from the perspective of the 'to' node
        if let Some(in_list) = self.edges_in.get_mut(&to_node) {
            in_list.push(edge);
        }
    }

    /// Generates unique k-mer anchors for the reference genome.
    /// `k` is the k-mer length (e.g., 31).
    pub fn build_reference_anchors(&mut self, sequence: &[u8], k: usize) -> Result<(), KmerError> {
        // A map to count the occurrences of every k-mer (both forward and RC)
        let mut kmer_counts: HashMap<u64, usize> = HashMap::new();

        // Pass 1: Count all k-mers in the combined forward and RC sets
        for kmer_slice in sequence.to_kmers(k)? {
            let fwd_encoded: u64 = kmer_slice.encode_to_u64()?;
            let rc_encoded:  u64 = kmer_slice.reverse_complement().encode_to_u64()?;

            *kmer_counts.entry(fwd_encoded).or_insert(0) += 1;
            *kmer_counts.entry(rc_encoded).or_insert(0) += 1;
        }

        // Pass 2: Identify forward k-mers that appear exactly once globally
        for (pos, kmer_slice) in sequence.to_kmers(k)?.enumerate() {
            let fwd_encoded: u64 = kmer_slice.encode_to_u64()?;

            // If the count is exactly 1, it means it appears once on the forward strand,
            // never on the reverse strand, and is not a palindrome (which would count as 2).
            if kmer_counts.get(&fwd_encoded) == Some(&1) {
                let anchor: ReferenceAnchor = ReferenceAnchor {
                    kmer: kmer_slice.to_vec(),
                    target_node_id: pos, // Node IDs match positions in the initial backbone
                };
                
                // Store using the 2-bit encoded u64 as the key for fast lookups later
                self.anchors.insert(fwd_encoded, anchor);
            }
        }

        Ok(())
    }

    /// Maps a query sequence to the reference graph using unique anchors.
    /// Returns the oriented query sequence and a sorted list of unique anchor matches.
    pub fn map_query_anchors(
        &self, 
        query_sequence: &[u8], 
        orientor: &StrandOrientor, 
        k: usize
    ) -> Result<(Vec<u8>, Vec<AnchorMatch>), KmerError> {
        
        // 1. Orient the query sequence using your StrandOrientor
        let oriented_query: Vec<u8> = orientor.orient(query_sequence)?;

        // 2. Scan the oriented query for reference anchors and track their query positions
        // We map the encoded k-mer to a Vec of its positions in the query
        let mut query_anchor_hits: HashMap<u64, Vec<usize>> = HashMap::new();

        for (pos, kmer_slice) in oriented_query.to_kmers(k)?.enumerate() {
            let encoded_kmer: u64 = kmer_slice.encode_to_u64()?;
            
            // If this k-mer is a known unique reference anchor, record its query position
            if self.anchors.contains_key(&encoded_kmer) {
                query_anchor_hits.entry(encoded_kmer).or_insert_with(Vec::new).push(pos);
            }
        }

        // 3. Filter out repeats in the query and build the final match list
        let mut unique_matches: Vec<AnchorMatch> = Vec::new();

        for (encoded_kmer, query_positions) in query_anchor_hits {
            // Strictly enforce that the anchor appears exactly once in the query genome
            if query_positions.len() == 1 {
                if let Some(ref_anchor) = self.anchors.get(&encoded_kmer) {
                    unique_matches.push(AnchorMatch {
                        query_position: query_positions[0],
                        reference_node_id: ref_anchor.target_node_id,
                    });
                }
            }
        }

        // 4. Sort the matches chronologically by their position in the query genome
        // This is crucial for comparing the lengths between anchors in the next step
        unique_matches.sort_by_key(|m| m.query_position);

        Ok((oriented_query, unique_matches))
    }

   /// Identifies the unaligned gaps between anchor matches.
    pub fn get_unaligned_segments(
        &self,
        matches: &[AnchorMatch],
        k: usize,
        query_len: usize,
    ) -> Vec<UnalignedSegment> {
        let mut segments: Vec<UnalignedSegment> = Vec::new();
        if matches.is_empty() { return segments; }

        let first_match: &AnchorMatch = &matches[0];
        let last_match:  &AnchorMatch = &matches[matches.len() - 1];

        // 1. 5' Flank
        if first_match.query_position > 0 && !self.is_circular {
            segments.push(UnalignedSegment {
                query_start: 0,
                query_end: first_match.query_position,
                segment_type: SegmentType::ExtendBackward { 
                    ref_node_end: first_match.reference_node_id 
                },
            });
        }

        // 2. Bounded Segments
        for window in matches.windows(2) {
            let m1: &AnchorMatch = &window[0];
            let m2: &AnchorMatch = &window[1];

            let q_start: usize = m1.query_position + k;
            let q_end:   usize = m2.query_position;

            if q_start < q_end || (m1.reference_node_id + k) != m2.reference_node_id {
                segments.push(UnalignedSegment {
                    query_start: q_start,
                    query_end: q_end,
                    segment_type: SegmentType::Bounded {
                        ref_node_start: m1.reference_node_id + k,
                        ref_node_end: m2.reference_node_id,
                    },
                });
            }
        }

        // 3. 3' Flank
        let last_anchor_end: usize = last_match.query_position + k;
        if last_anchor_end < query_len && !self.is_circular {
            segments.push(UnalignedSegment {
                query_start: last_anchor_end,
                query_end: query_len,
                segment_type: SegmentType::ExtendForward { 
                    ref_node_start: last_match.reference_node_id + k 
                },
            });
        }

        segments
    }

    /// Extracts a topologically sorted list of node IDs between a start and end node.
    /// This is required to run Partial Order Alignment on a specific gap.
    pub fn get_topological_subgraph(&self, start_node: usize, end_node: usize) -> Vec<usize> {
        let mut in_degrees: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
        let mut queue: std::collections::VecDeque<usize> = std::collections::VecDeque::new();
        let mut sorted_nodes: Vec<usize> = Vec::new();

        // 1. Discover all nodes in this subgraph and count their local in-degrees
        // We use a basic BFS starting from `start_node`.
        let mut to_visit: Vec<usize> = vec![start_node];
        let mut visited: std::collections::HashSet<usize> = std::collections::HashSet::new();
        visited.insert(start_node);

        while let Some(current) = to_visit.pop() {
            // Stop traversing down this path if we hit the boundary
            if current == end_node {
                continue; 
            }

            if let Some(edges) = self.edges_out.get(&current) {
                for edge in edges {
                    let next_node = edge.to_node;
                    
                    // Increment the in-degree for this local subgraph
                    *in_degrees.entry(next_node).or_insert(0) += 1;

                    if !visited.contains(&next_node) {
                        visited.insert(next_node);
                        to_visit.push(next_node);
                    }
                }
            }
        }

        // 2. Kahn's Algorithm for Topological Sorting
        // Start with the start_node, which has an artificial local in-degree of 0
        queue.push_back(start_node);

        while let Some(current) = queue.pop_front() {
            sorted_nodes.push(current);

            if current == end_node {
                continue; // Don't traverse past the boundary
            }

            if let Some(edges) = self.edges_out.get(&current) {
                for edge in edges {
                    let next_node = edge.to_node;
                    
                    if let Some(degree) = in_degrees.get_mut(&next_node) {
                        *degree -= 1;
                        if *degree == 0 {
                            queue.push_back(next_node);
                        }
                    }
                }
            }
        }

        sorted_nodes
    }
}

// ============================================================================
// Alignment graph - Align segments
// ============================================================================

impl AlignmentGraph {
    /// Weaves the results of a dynamic programming alignment into the graph.
    pub fn weave_alignment(
        &mut self,
        ops:          &[EditOp],
        ref_path:     &[usize],
        anchor_start: Option<usize>,
        anchor_end:   Option<usize>,
    ) {
        let mut ref_idx: usize = 0;
        let mut last_query_node: Option<usize> = anchor_start;

        for op in ops {
            match op {
                EditOp::Match => {
                    let current_ref_node: usize = ref_path[ref_idx];
                    
                    if let Some(prev) = last_query_node {
                        self.add_edge(prev, Strand::Forward, current_ref_node, Strand::Forward, EdgeType::Match);
                    }
                    
                    last_query_node = Some(current_ref_node);
                    ref_idx += 1;
                }
                EditOp::Substitution(base) => {
                    let new_node: usize = self.add_position_node(0, *base);
                    if let Some(prev) = last_query_node {
                        self.add_edge(prev, Strand::Forward, new_node, Strand::Forward, EdgeType::Mismatch);
                    }
                    last_query_node = Some(new_node);
                    ref_idx += 1;
                }
                EditOp::Insertion(base) => {
                    let new_node: usize = self.add_position_node(0, *base);
                    if let Some(prev) = last_query_node {
                        self.add_edge(prev, Strand::Forward, new_node, Strand::Forward, EdgeType::Insertion);
                    }
                    last_query_node = Some(new_node);
                }
                EditOp::Deletion => {
                    ref_idx += 1;
                }
            }
        }

        if let (Some(last), Some(end_node)) = (last_query_node, anchor_end) {
            self.add_edge(last, Strand::Forward, end_node, Strand::Forward, EdgeType::Match);
        }
    }

    /// Prepares the sorted nodes for the POA algorithm by mapping their bases and predecessors.
    fn build_poa_inputs(&self, sorted_nodes: &[usize]) -> (Vec<u8>, Vec<Vec<usize>>) {
        let mut node_bases: Vec<u8> = Vec::with_capacity(sorted_nodes.len());
        let mut predecessors: Vec<Vec<usize>> = vec![Vec::new(); sorted_nodes.len()];
        
        // Map real node IDs to their topological array index (j)
        let mut node_to_idx = std::collections::HashMap::new();
        for (i, &node_id) in sorted_nodes.iter().enumerate() {
            node_bases.push(self.positions.get(&node_id).expect("Node missing").nucleotide);
            node_to_idx.insert(node_id, i);
        }
        
        // Build the predecessor array
        for (i, &node_id) in sorted_nodes.iter().enumerate() {
            if let Some(in_edges) = self.edges_in.get(&node_id) {
                for edge in in_edges {
                    // Only include predecessors that are actually inside this subgraph
                    if let Some(&p_idx) = node_to_idx.get(&edge.from_node) {
                        predecessors[i].push(p_idx);
                    }
                }
            }
        }
        
        (node_bases, predecessors)
    }

    /// Weaves a POA alignment result into the graph.
    pub fn weave_poa_alignment(
        &mut self,
        ops: &[EditOp],
        poa_path: &[usize],     // Indices corresponding to traversed nodes in sorted_nodes
        sorted_nodes: &[usize], // The actual graph node IDs
        anchor_start: Option<usize>,
        anchor_end: Option<usize>,
    ) {
        let mut path_idx = 0;
        let mut last_query_node = anchor_start;

        for op in ops {
            match op {
                EditOp::Match => {
                    let current_ref_node = sorted_nodes[poa_path[path_idx]];
                    if let Some(prev) = last_query_node {
                        self.add_edge(prev, Strand::Forward, current_ref_node, Strand::Forward, EdgeType::Match);
                    }
                    last_query_node = Some(current_ref_node);
                    path_idx += 1;
                }
                EditOp::Substitution(base) => {
                    let new_node = self.add_position_node(0, *base);
                    if let Some(prev) = last_query_node {
                        self.add_edge(prev, Strand::Forward, new_node, Strand::Forward, EdgeType::Mismatch);
                    }
                    last_query_node = Some(new_node);
                    path_idx += 1; // It replaces a graph node, so we advance the graph cursor
                }
                EditOp::Insertion(base) => {
                    let new_node = self.add_position_node(0, *base);
                    if let Some(prev) = last_query_node {
                        self.add_edge(prev, Strand::Forward, new_node, Strand::Forward, EdgeType::Insertion);
                    }
                    last_query_node = Some(new_node);
                    // Insertion does not consume a graph node
                }
                EditOp::Deletion => {
                    path_idx += 1; // Advance past the deleted graph node silently
                }
            }
        }

        if let (Some(last), Some(end_node)) = (last_query_node, anchor_end) {
            self.add_edge(last, Strand::Forward, end_node, Strand::Forward, EdgeType::Match);
        }
    }
}

// ============================================================================
// Alignment graph - Build
// ============================================================================

impl AlignmentGraph {
    /// Full pipeline to initialize a reference graph and map a query genome onto it.
    pub fn build_pangenome(
        reference_seq: &[u8],
        query_sequences: &[&[u8]], // A list of query genomes (e.g., from a parsed Multi-FASTA)
        k: usize,
        is_circular: bool,
        export_path: &str,
    ) -> Result<AlignmentGraph, Box<dyn std::error::Error>> {
        
        println!("Initializing graph and adding reference backbone...");
        let mut graph: AlignmentGraph = AlignmentGraph::new(is_circular)?;
        graph.add_reference(reference_seq);

        println!("Building reference anchors (k={})...", k);
        graph.build_reference_anchors(reference_seq, k)?;

        let orientor: StrandOrientor = StrandOrientor::new(reference_seq, 15)?;
        let scores: AlignmentScores  = AlignmentScores::default();

        for (i, query_seq) in query_sequences.iter().enumerate() {
            println!("Aligning Genome #{}...", i + 2); // Reference is #1
            graph.add_genome(query_seq, k, &orientor, &scores)?;
        }

        println!("Exporting graph to {}...", export_path);
        graph.export_to_graphml(export_path)?;

        println!("Pangenome construction complete!");
        Ok(graph)
    }
}

// ============================================================================
// Alignment graph - IO
// ============================================================================

impl AlignmentGraph {
    /// Exports the graph to Cytoscape-friendly GraphML
    pub fn export_to_graphml(&self, path: &str) -> io::Result<()> {
        let file: File = File::create(path)?;
        let mut writer: BufWriter<File> = BufWriter::new(file);

        writeln!(writer, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>")?;
        writeln!(writer, "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\">")?;
        
        // Define Node attributes
        writeln!(writer, "  <key id=\"v_type\" for=\"node\" attr.name=\"type\" attr.type=\"string\"/>")?;
        writeln!(writer, "  <key id=\"v_label\" for=\"node\" attr.name=\"label\" attr.type=\"string\"/>")?;
        writeln!(writer, "  <key id=\"v_pos\" for=\"node\" attr.name=\"ref_pos\" attr.type=\"int\"/>")?;
        // (Note: You can add 'frequency' back here later when we track coverage!)

        // Define Edge attributes
        writeln!(writer, "  <key id=\"e_type\" for=\"edge\" attr.name=\"type\" attr.type=\"string\"/>")?;
        writeln!(writer, "  <key id=\"e_from_strand\" for=\"edge\" attr.name=\"from_strand\" attr.type=\"string\"/>")?;
        writeln!(writer, "  <key id=\"e_to_strand\" for=\"edge\" attr.name=\"to_strand\" attr.type=\"string\"/>")?;

        writeln!(writer, "  <graph id=\"alignment_graph\" edgedefault=\"directed\">")?;

        // 1. Write Nodes
        // We iterate over the values in our positions HashMap
        for node in self.positions.values() {
            writeln!(writer, "    <node id=\"n{}\">", node.id)?;
            writeln!(writer, "      <data key=\"v_type\">Nucleotide</data>")?;
            writeln!(writer, "      <data key=\"v_label\">{}</data>", node.nucleotide as char)?;
            writeln!(writer, "      <data key=\"v_pos\">{}</data>", node.position)?;
            writeln!(writer, "    </node>")?;
        }

        // 2. Write Edges
        // We flatten our edges_out HashMap to get all outgoing edges
        for edge in self.edges_out.values().flatten() {
            let e_type_str = match edge.edge_type {
                EdgeType::Backbone  => "Backbone",
                EdgeType::Match     => "Match",
                EdgeType::Mismatch  => "Mismatch",
                EdgeType::Insertion => "Insertion",
                EdgeType::Deletion  => "Deletion", // If explicitly tracked
            };

            let from_strand_str = match edge.from_strand {
                Strand::Forward => "Forward",
                Strand::Reverse => "Reverse",
            };

            let to_strand_str = match edge.to_strand {
                Strand::Forward => "Forward",
                Strand::Reverse => "Reverse",
            };

            writeln!(writer, "    <edge source=\"n{}\" target=\"n{}\">", edge.from_node, edge.to_node)?;
            writeln!(writer, "      <data key=\"e_type\">{}</data>", e_type_str)?;
            writeln!(writer, "      <data key=\"e_from_strand\">{}</data>", from_strand_str)?;
            writeln!(writer, "      <data key=\"e_to_strand\">{}</data>", to_strand_str)?;
            writeln!(writer, "    </edge>")?;
        }

        writeln!(writer, "  </graph>\n</graphml>")?;
        Ok(())
    }
}


// ============================================================================
// Graph and alignment components
// ============================================================================

/// Represents the reading direction of a node.
#[derive(Debug, Clone)]
pub enum Strand {
    Forward, // Reading 5' to 3'
    Reverse, // Reading 3' to 5' (requires reverse complementing the nucleotide)
}

impl Strand {
    /// Flips the strand direction.
    pub fn flip(&self) -> Self {
        match self {
            Strand::Forward => Strand::Reverse,
            Strand::Reverse => Strand::Forward,
        }
    }
}

/// Represents a 1-to-1 anchor match between the query and the reference.
#[derive(Debug, Clone)]
pub struct AnchorMatch {
    pub query_position: usize,
    pub reference_node_id: usize,
}

#[derive(Debug, Clone)]
pub struct ReferencePosition {
    pub id: usize,
    pub position: usize,
    pub nucleotide: u8,
}

#[derive(Debug, Clone)]
pub struct ReferenceAnchor {
    pub kmer: Vec<u8>,
    pub target_node_id: usize,
}

#[derive(Debug, Clone)]
pub enum SegmentType {
    /// Bounded by two anchors. Aligned normally (left-to-right).
    Bounded { 
        ref_node_start: usize, 
        ref_node_end: usize 
    },
    
    /// 3' Flank (Linear only). Anchored at the left, extends right.
    ExtendForward { 
        ref_node_start: usize 
    },
    
    /// 5' Flank (Linear only). Anchored at the right, extends left.
    /// Sequences must be reversed (not complemented!) before alignment.
    ExtendBackward { 
        ref_node_end: usize 
    },
}

/// Represents a slice of the query genome that falls between two anchors
/// and requires more granular alignment against the reference.
#[derive(Debug, Clone)]
pub struct UnalignedSegment {
    pub query_start: usize,
    pub query_end: usize,
    pub segment_type: SegmentType,
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    /// Helper function to initialize a basic graph with a linear backbone
    fn setup_mock_graph(reference: &[u8]) -> AlignmentGraph {
        let mut graph: AlignmentGraph = AlignmentGraph {
            positions: std::collections::HashMap::new(),
            edges_out: std::collections::HashMap::new(),
            edges_in: std::collections::HashMap::new(),
            anchors: std::collections::HashMap::new(),
            next_node_id: 0,
            is_circular: false,
        };
        
        // Build backbone manually for testing
        for (i, &base) in reference.iter().enumerate() {
            let node_id = graph.next_node_id;
            graph.next_node_id += 1;
            
            graph.positions.insert(node_id, ReferencePosition {
                id: node_id,
                nucleotide: base,
                position: i,
            });
            
            if i > 0 {
                // Link previous node to this one
                graph.add_edge(
                    node_id - 1, Strand::Forward, 
                    node_id, Strand::Forward, 
                    EdgeType::Backbone
                );
            }
        }
        graph
    }

    #[test]
    fn test_snp_bubble_creation() {
        // Ref: A T G C A
        // Qry: A T A C A (SNP at pos 2: G -> A)
        let ref_seq: &[u8; 5] = b"ATGCA";
        let query_seq: &[u8; 5] = b"ATACA";
        let mut graph: AlignmentGraph = setup_mock_graph(ref_seq);

        let scores: AlignmentScores = AlignmentScores::default();
        let ops: Vec<EditOp> = align_seqs(ref_seq, query_seq, &scores);

        assert_eq!(ops, vec![
            EditOp::Match,
            EditOp::Match,
            EditOp::Substitution(b'A'),
            EditOp::Match,
            EditOp::Match,
        ]);

        // Weave the alignment. 
        // We branch off a hypothetical anchor at index 0, and merge at index 4
        let ref_path: Vec<usize> = (0..ref_seq.len()).collect();
        graph.weave_alignment(&ops, &ref_path, Some(0), Some(4));

        // Export to visually inspect in Cytoscape
        let export_path: &str = "test_snp_bubble.graphml";
        graph.export_to_graphml(export_path).unwrap();

        // Check file was created
        assert!(fs::metadata(export_path).is_ok());
        
        // Clean up test file
        let _ = fs::remove_file(export_path);
    }

    #[test]
    fn test_insertion_bubble() {
        // Ref: A T C A
        // Qry: A T G C A (Insertion of G)
        let ref_seq: &[u8; 4] = b"ATCA";
        let query_seq: &[u8; 5] = b"ATGCA";
        let mut graph: AlignmentGraph = setup_mock_graph(ref_seq);

        let scores: AlignmentScores = AlignmentScores::default();
        let ops: Vec<EditOp> = align_seqs(ref_seq, query_seq, &scores);

        assert_eq!(ops, vec![
            EditOp::Match,
            EditOp::Match,
            EditOp::Insertion(b'G'),
            EditOp::Match,
            EditOp::Match,
        ]);

        let ref_path: Vec<usize> = (0..ref_seq.len()).collect();
        graph.weave_alignment(&ops, &ref_path, Some(0), Some(3));

        let export_path: &str = "test_insertion_bubble.graphml";
        graph.export_to_graphml(export_path).unwrap();
        let _ = fs::remove_file(export_path);
    }

    #[test]
    fn test_deletion_shortcut() {
        // Ref: A T G C A
        // Qry: A T - C A (Deletion of G)
        let ref_seq: &[u8; 5] = b"ATGCA";
        let query_seq: &[u8; 4] = b"ATCA";
        let mut graph: AlignmentGraph = setup_mock_graph(ref_seq);

        let scores: AlignmentScores = AlignmentScores::default();
        let ops: Vec<EditOp> = align_seqs(ref_seq, query_seq, &scores);

        assert_eq!(ops, vec![
            EditOp::Match,
            EditOp::Match,
            EditOp::Deletion,
            EditOp::Match,
            EditOp::Match,
        ]);

        let ref_path: Vec<usize> = (0..ref_seq.len()).collect();
        graph.weave_alignment(&ops, &ref_path, Some(0), Some(4));

        let export_path: &str = "test_deletion_shortcut.graphml";
        graph.export_to_graphml(export_path).unwrap();
        let _ = fs::remove_file(export_path);
    }

    #[test]
    fn test_topological_subgraph() {
        let mut graph: AlignmentGraph = setup_mock_graph(b"A"); // Dummy backbone
        
        // We will manually construct a bubble:
        //      /-> 2 (C) -\
        // 1 (A)            -> 4 (T)
        //      \-> 3 (G) -/
        
        // Nodes 1 through 4
        let n1: usize = graph.add_position_node(1, b'A');
        let n2: usize = graph.add_position_node(2, b'C');
        let n3: usize = graph.add_position_node(3, b'G');
        let n4: usize = graph.add_position_node(4, b'T');
        
        // Edges
        graph.add_edge(n1, Strand::Forward, n2, Strand::Forward, EdgeType::Match);
        graph.add_edge(n1, Strand::Forward, n3, Strand::Forward, EdgeType::Mismatch);
        graph.add_edge(n2, Strand::Forward, n4, Strand::Forward, EdgeType::Match);
        graph.add_edge(n3, Strand::Forward, n4, Strand::Forward, EdgeType::Match);

        let sorted: Vec<usize> = graph.get_topological_subgraph(n1, n4);

        // 1 must be first, 4 must be last. 2 and 3 can be in any order.
        assert_eq!(sorted.first(), Some(&n1));
        assert_eq!(sorted.last(), Some(&n4));
        assert_eq!(sorted.len(), 4);
        assert!(sorted.contains(&n2));
        assert!(sorted.contains(&n3));
    }
}