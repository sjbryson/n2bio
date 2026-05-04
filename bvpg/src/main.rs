//! n2bio/bvpg/src/main.rs

use clap::Parser;
use serde::Serialize;
use std::process;
use n2core::kgraph::PanGenomeGraph;


#[derive(Parser, Debug)]
#[command(name = "PanGenomeBuilder")]
#[command(author, version, about = "Builds a directed Pan-Genome Graph from FASTA sequences and exports to GraphML for Cytoscape.")]
struct Args {
    /// Path to the reference FASTA file (backbone)
    #[arg(long, required = true)]
    reference: String,

    /// Path to the multi-FASTA assemblies file
    #[arg(long, required = true)]
    genomes: String,

    /// K-mer size (Must be > 0 and <= 32 for 2-bit encoding)
    #[arg(short = 'k', long = "kmer_size", default_value_t = 31)]
    kmer_size: usize,

    /// Output path for the GraphML file
    #[arg(long, required = true)]
    graphml: String,
}

/// Struct to hold graph stats
#[derive(Serialize)]
struct GraphStats<'a> {
    kmer_size: usize,
    total_nodes: usize,
    total_edges: usize,
    reference_file: &'a str,
    assemblies_file: &'a str,
    output_graphml: &'a str,
    status: &'static str,
}

fn main() {
    // 1. Parse command line arguments
    let args = Args::parse();

    // 2. Validate k-mer size
    if args.kmer_size == 0 || args.kmer_size > 32 {
        eprintln!("Error: k-mer size (-k) must be between 1 and 32.");
        process::exit(1);
    }

    // 3. Build the graph
    match PanGenomeGraph::from_fastas(&args.reference, &args.genomes, args.kmer_size) {
        Ok(graph) => {
            // 4. Export to GraphML
            if let Err(e) = graph.export_to_graphml(&args.graphml) {
                eprintln!("Failed to export GraphML file: {}", e);
                process::exit(1);
            }

            // 5. JSON formatted graph info to stdout
            let stats = GraphStats {
                kmer_size: args.kmer_size,
                total_nodes: graph.graph.node_count(),
                total_edges: graph.graph.edge_count(),
                reference_file: &args.reference,
                assemblies_file: &args.genomes,
                output_graphml: &args.graphml,
                status: "success",
            };

            match serde_json::to_string_pretty(&stats) {
                Ok(json_output) => println!("{}", json_output),
                Err(e) => {
                    eprintln!("Failed to serialize JSON stats: {}", e);
                    process::exit(1);
                }
            }
        }
        Err(e) => {
            eprintln!("Failed to build graph due to an IO error: {}", e);
            process::exit(1);
        }
    }
}