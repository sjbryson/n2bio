//! n2bio/bvgg/src/main.rs

use clap::Parser;
use std::process;
use n2core::alignmentgraph::AlignmentGraph;


#[derive(Parser, Debug)]
#[command(name = "PanGenomeBuilder")]
#[command(author, version, about = "Builds a directed AlignmentGraph from FASTA sequences and exports to GraphML for Cytoscape.")]
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

    /// Reference is circular (bool)
    #[arg(long, default_value_t = false)]
    circular: bool,

    /// Output path for the GraphML file
    #[arg(long, required = true)]
    graphml: String,
}

fn main() {
    // 1. Parse command line arguments
    let args: Args = Args::parse();

    // 2. Validate k-mer size
    if args.kmer_size == 0 || args.kmer_size > 32 {
        eprintln!("Error: k-mer size (-k) must be between 1 and 32.");
        process::exit(1);
    }

    // 3. Build the graph
    match AlignmentGraph::from_fastas(&args.reference, &args.genomes, args.kmer_size, args.circular, &args.graphml) {
        Ok(graph) => {
            // 4. Export to GraphML
            if let Err(e) = graph.export_to_graphml(&args.graphml) {
                eprintln!("Failed to export GraphML file: {}", e);
                process::exit(1);
            }
        }
        Err(e) => {
            eprintln!("Failed to build graph due to an IO error: {}", e);
            process::exit(1);
        }
    }
}