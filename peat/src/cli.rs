//! n2bio/peat/src/cli.rs
//! 

#![allow(unused)]

use clap::{ Args, Parser, Subcommand, ValueEnum };
use std::path::PathBuf;

// ============================================================================
// Subcommands
// ============================================================================

#[derive(Parser)]
#[command(name = "peat", version = "1.0", about = "Paired-End Alignment Tools")]
pub(crate) struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub(crate) enum Commands {
    /// Parse SAM records from stdin and filter to create a filtered paired-end fastq.gz library (r1.fq.gz & r2.fq.gz)
    Filter(FilterArgs),
    /// Parse SAM records from stdin and calculate coverage for each reference in the sam/bam header
    Coverage(CoverageArgs),
    /// Read a name sorted bam file and generate an interactive report
    BamRep(BamRepArgs),
    /// ToDo: Parse SAM records from stdin and bin read pairs for each target
    BinReads(BinReadsArgs),
}

// ============================================================================
// Threshold Mode
// ============================================================================

#[derive(ValueEnum, Clone, Debug)]
pub(crate) enum ThresholdMode {
    /// Low Pass: process SAM records that are below defined thresholds
    #[value(name = "lowpass")]
    LowPass,
    
    /// High Pass: process SAM records that are above defined thresholds
    #[value(name = "highpass")]
    HighPass,
}

// ============================================================================
// Threshold Metrics
// ============================================================================

#[derive(Args, Clone)]
pub(crate) struct ThresholdMetrics {
    /// Optional: Alignment Score - sam.get_int_tag("AS")
    #[arg(long = "align_score", alias = "AS")]
    pub align_score: Option<i32>,
    
    /// Optional: Alignment Lenth - sam.calculate_alignment_length()
    #[arg(long = "align_length", alias = "AL")]
    pub align_length: Option<u32>,

    /// Optional: Per base alignment score (AS/AL = avg. align_score per covered base) - sam.calculate_as_al()
    #[arg(long = "base_score", alias = "BS")]
    pub base_score: Option<f32>,

    /// Optional: Alignment Proportion - sam.calculate_alignment_proportion()
    #[arg(long = "align_prop", alias = "AP")]
    pub align_prop: Option<f32>,
    
    /// Optional: Alignment Percent Identity - sam.calculate_alignment_accuracy()
    #[arg(long = "align_ident", alias = "AI")]
    pub align_ident: Option<f32>,

    /// Optional: Max MAPQ score - sam.mapq()
    #[arg(long = "mapq", alias = "MQ")]
    pub mapq: Option<u32>,
}

// ============================================================================
// Filter Args
// ============================================================================

#[derive(Args, Clone)]
pub(crate) struct FilterArgs {
    
     /// Number of worker threads for parsing and pairing
    #[arg(short = 't', long, default_value_t = 4)]
    pub threads: usize,

    /// Number of shards for the mate pairing hash map (recommend 4-8x threads, default = 32)
    #[arg(long, default_value_t = 32)]
    pub shards: usize,

    /// Prefix for output files (e.g. 'out' -> out.r1.fq.gz, out.r2.fq.gz)
    #[arg(short = 'p', long = "prefix", required = true)]
    pub prefix: String,
    
    // make optional if --stdout interleaved or --interleaved (pe vs lr) options

    /// Name of the run/sample for the JSON report -> creates {report}.json
    #[arg(short = 'r', long = "report", required = true)]                              // ====== Update output logic ===========
    pub report: String,

   /// How abundance values should be mathematically interpreted
    #[arg(long = "filter_mode", value_enum)]
    pub filter_mode: ThresholdMode,

    #[command(flatten)]
    pub thresholds: ThresholdMetrics
    
}

// ============================================================================
// Coverage Args
// ============================================================================

#[derive(Args, Clone)]
pub(crate) struct CoverageArgs {
    
    /// Number of worker threads for parsing
    #[arg(short = 't', long, default_value_t = 4)]
    pub threads: usize,

    /// Name of the run/sample for the JSON report -> creates {report}.json
    #[arg(short = 'r', long, required = true)]
    pub report: String,

     /// Optional path to an SQLite taxonomy database (see vref2db)
    #[arg(long)]
    pub db: Option<String>,

    #[command(flatten)]
    pub thresholds: ThresholdMetrics
}

// ============================================================================
// Compose Args
// ============================================================================

#[derive(Args)]
pub(crate) struct BamRepArgs {
    
    /// Input name-sorted BAM file
    #[arg(short = 'b', long)]
    pub bam: PathBuf,

    /// Report file prefix - creates {report}.json and optional {report}.html
    #[arg(short = 'r', long)]
    pub report: PathBuf,

    /// Generate html plots
    #[arg(long)]
    pub html: bool,

    /// Minimum MAPQ score for insert size calculation
    #[arg(short = 'q', long, default_value_t = 40)]
    pub min_mapq: usize,

    /// Max insert size to use for summary stats calculation
    #[arg(short = 'i', long, default_value_t = 1000)]
    pub max_ins: usize,

    /// Max read length to use
    #[arg(short = 'l', long, default_value_t = 150)]
    pub max_len: usize,
}

// ============================================================================
// BinReads Args
// ============================================================================

#[derive(Parser, Debug, Clone)]
pub(crate) struct BinReadsArgs {
    
    /// Path to an input BAM file to evaluate
    #[arg(short = 'b', long)]
    pub bam: String,

    /// Path to a TSV mapping file: referenc_ id --> bin_id
    #[arg(short = 'r', long = "reference-map")]
    pub reference_map: String,
}

