//! n2bio/peat/src/cli.rs
//! 

#![allow(unused)]

use clap::{ Args, Parser, Subcommand };

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
    /// 
    Coverage(CoverageArgs),
    /// 
    BamRep(BamRepArgs),
    //
    BinReads(BinReadsArgs),
}

// ============================================================================
// Threshold Mode
// ============================================================================

#[derive(clap::ValueEnum, Clone, Debug)]
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

#[derive(Args)]
pub(crate) struct CoverageArgs {
    
    /// Path to the JSON model report -> created by pfqsim model
    #[arg(short = 'm', long)]
    pub model: String,

    /// Path to the genome fasta file to generate reads from
    #[arg(short = 'f', long)]
    pub fasta: String,

    /// Boolean value: circularize genome before generating reads
    #[arg(short = 'c', long, default_value_t = false)]
    pub circular: bool,

    /// Float value for random substitution rate to apply to simulated reads (range: 0.0 - 1.0)
    #[arg(short = 's', long, default_value_t = 0.0)]
    pub sub_rate: f64,

    /// Float value for random insertion and deletion rate to apply to simulated reads (range: 0.0 - 1.0)
    #[arg(short = 'i', long, default_value_t = 0.0)]
    pub indel_rate: f64,

    /// Integer value for number of paired reads to create (1 = 1 R1.fq.gz + 1 R2.fq.gz)
    #[arg(short = 'n', long)]
    pub num_reads: usize,

    /// Read length to generate
    #[arg(short = 'l', long, default_value_t = 150)]
    pub read_length: usize,

    /// Boolean value: Vary read lengths based on model
    #[arg(long, default_value_t = false)]
    pub vary_lengths: bool,

    /// Prefix for output fastq.gz files (e.g. {prefix}.r1.fq.gz)
    /// and for read identifiers (e.g. @{prefix}:{keyword}:Accession::Read Num)
    #[arg(short = 'p', long)]
    pub prefix: String,

    /// Additional keyword to add to read identifiers
    /// for use in query-target mapping (e.g. @{prefix}:{keyword}:Accession::Read Num)
    #[arg(short = 'k', long)]
    pub keyword: String,

    /// Number of worker threads
    #[arg(short = 't', long)]
    pub threads: usize,

    #[arg(skip)]
    pub append_path: Option<String>,

    #[arg(skip)]
    pub append_mode: bool,

}

// ============================================================================
// Compose Args
// ============================================================================

#[derive(Args)]
pub(crate) struct BamRepArgs {
    
    /// Path to a TSV config file 
    #[arg(short = 'c', long)]
    pub config: String,

    /// Path to the JSON model report -> created by pfqsim model
    #[arg(short = 'm', long)]
    pub model: String,

    /// Prefix for the manifest tsv and both simulated reads (R1 & R2) files 
    #[arg(short = 'p', long)]
    pub prefix: String,

    /// Integer value for number of paired reads to create (1 = 1 R1.fq.gz + 1 R2.fq.gz)
    #[arg(short = 'n', long)]
    pub total_reads: usize,

    /// Read length to generate
    #[arg(short = 'l', long, default_value_t = 150)]
    pub read_length: usize,

    /// Boolean value: Vary read lengths based on model
    #[arg(long, default_value_t = false)]
    pub vary_lengths: bool,

    /// Number of worker threads
    #[arg(short = 't', long)]
    pub threads: usize,
}

// ============================================================================
// BinReads Args
// ============================================================================

#[derive(Parser, Debug, Clone)]
pub(crate) struct BinReadsArgs {
    
    /// Path to a TSV config file 
    #[arg(short = 'c', long)]
    pub config: String,

    /// Path to an input BAM file to evaluate
    #[arg(short = 'b', long)]
    pub bam: String,

    /// Path to a TSV mapping file: reference id --> mapping-mode(id, keyword, or accession)
    #[arg(short = 'r', long = "reference-map")]
    pub reference_map: String,
}

