//! n2bio/pfqsim/src/cli.rs
//! 

use clap::{ Args, Parser, Subcommand };
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "pfqsim", version = "1.0", about = "Fast metagenomic read simulator")]
pub(crate) struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub(crate) enum Commands {
    /// Build insert size and Q-score distributions from a BAM file
    Model(ModelArgs),
    /// Generate a simulated paired-read library from a reference FASTA
    Generate(GenerateArgs),
    /// Compose a final metagenomic library based on an abundance config
    Compose(ComposeArgs),
    /// Analyze alignments from stdin sam or a bam file
    Analyze(AnalyzeArgs),
}

#[derive(Args)]
pub(crate) struct ModelArgs {
    
    /// Name for the BAM file for modeling insert size and Qscore distributions
    #[arg(short = 'b', long)]
    pub bam: PathBuf,

    /// Name for the JSON model report -> creates {model}.json
    #[arg(short = 'm', long)]
    pub model: PathBuf,

    /// Read length to model (default = 150)
    #[arg(short = 'l', long, default_value_t = 150)]
    pub read_length: usize,

    /// Min MAPQ score for filtering alignments for insert size distribution
    #[arg(short = 'q', long, default_value_t = 40)]
    pub mapq: usize,

    /// Max insert size to use for insert size distribution
    #[arg(short = 'i', long, default_value_t = 1000)]
    pub max_ins: usize,
}

#[derive(Args)]
pub(crate) struct GenerateArgs {
    
    /// Path to the JSON model report -> created by pfqsim --model
    #[arg(short = 'm', long)]
    pub model: PathBuf,

    /// Path to the fasta file to generate reads from
    #[arg(short = 'f', long)]
    pub fasta: PathBuf,

    /// Boolean value: circular genome
    #[arg(short = 'c', long, default_value_t = false)]
    pub circular: bool,

    /// Float value for random substitution rate to apply to simulated reads
    #[arg(short = 's', long)]
    pub sub_rate: f64,

    /// Float value for random insertion and deletion rate to apply to simulated reads
    #[arg(short = 'i', long)]
    pub indel_rate: f64,

    /// Integer value for number of paired reads to create (1 = 1 R1.fq.gz + 1 R2.fq.gz)
    #[arg(short = 'n', long)]
    pub num_reads: usize,

    /// Read length to model (default = 150)
    #[arg(short = 'l', long, default_value_t = 150)]
    pub read_length: usize,

    /// Prefix for output fastq.gz files (e.g. {prefix}.r1.fq.gz)
    /// and for read identifiers (e.g. @{prefix}:Accession::Read Num)
    #[arg(short = 'p', long)]
    pub prefix: String,

    /// Number of worker threads
    #[arg(short = 't', long)]
    pub threads: usize,
}


#[derive(clap::ValueEnum, Clone, Copy, Debug, PartialEq)]
pub(crate) enum AbundanceMode {
    /// Calculate abundance as fraction of total reads
    #[value(name = "reads")]
    ReadFraction,
    
    /// Calculate abundance as fraction of total genome copies
    #[value(name = "copies")]
    CopyFraction,
}

#[derive(Args)]
pub(crate) struct ComposeArgs {
    
    /// Path to a TSV config file 
    #[arg(short = 'c', long)]
    pub config: PathBuf,

    #[arg(short = 'p', long)]
    pub prefix: PathBuf,

    /// Integer value for number of paired reads to create (1 = 1 R1.fq.gz + 1 R2.fq.gz)
    #[arg(short = 'n', long)]
    pub total_reads: usize,

    /// Number of worker threads
    #[arg(short = 't', long)]
    pub threads: usize,

    /// How abundance values should be mathematically interpreted
    #[arg(long = "abundance-mode", value_enum, default_value_t = AbundanceMode::ReadFraction)]
    pub abundance_mode: AbundanceMode,
}

#[derive(Args)]
pub(crate) struct AnalyzeArgs {
    
    /// Path to a TSV config file 
    #[arg(short = 'c', long)]
    pub config: PathBuf,

    /// Name of the BAM file for analyzing
    #[arg(short = 'b', long)]
    pub bam: PathBuf,

    /// Name for the analysis report
    #[arg(short = 'o', long)]
    pub output: PathBuf,

    // ToDo: option for stdin sam input
}