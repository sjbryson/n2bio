//! n2bio/pfqsim/src/cli.rs
//! 

use clap::{ Args, Parser, Subcommand };

#[derive(Parser)]
#[command(name = "pfqsim", version = "1.0", about = "Fast metagenomic read simulator")]
pub(crate) struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub(crate) enum Commands {
    /// Build insert size, read length, and Q-score distributions from a name sorted BAM file
    Model(ModelArgs),
    /// Generate a simulated paired-read library from a reference FASTA
    Generate(GenerateArgs),
    /// Compose a final metagenomic library based on an abundance config
    Compose(ComposeArgs),
    /// Analyze alignments from a name sorted BAM file for classification stats
    Analyze(AnalyzeArgs),
}

#[derive(Args)]
pub(crate) struct ModelArgs {
    
    /// Path to the BAM file for modeling insert size, read length, and Qscore distributions
    #[arg(short = 'b', long)]
    pub bam: String,

    /// Name for the JSON model report -> creates {model}.json
    #[arg(short = 'm', long)]
    pub model: String,

    /// Read length to model
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


#[derive(clap::ValueEnum, Clone, Debug)]
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

    /// How abundance values should be mathematically interpreted
    #[arg(long = "abundance-mode", value_enum, default_value_t = AbundanceMode::ReadFraction)]
    pub abundance_mode: AbundanceMode,
}

#[derive(Parser, Debug, Clone)]
pub(crate) struct AnalyzeArgs {
    
    /// Path to a TSV config file 
    #[arg(short = 'c', long)]
    pub config: String,

    /// Path to an input BAM file to evaluate
    #[arg(short = 'b', long)]
    pub bam: String,

    /// Path to a TSV mapping file: reference id --> mapping-mode(id, keyword, or accession)
    #[arg(short = 'r', long = "reference-map")]
    pub reference_map: String,

    /// Output path prefix for the generated HTML & json evaluation reports
    #[arg(short = 'o', long, default_value = "pfqsim_report")]
    pub output: String,

    /// Part of read identifier to use for reference sequence mapping
    #[arg(short = 'm', long = "mapping-mode")]
    pub mapping_mode: MappingMode,
}

/// Which part of the read identifier to map to database accessions for analysis
/// @{ReadId}:{ReadKeyword}:{ReadAccession}:ReadNumber
#[derive(clap::ValueEnum, Clone, Debug)]
pub(crate) enum MappingMode {
    /// Reference database mapping: Accession --> ReadId
    #[value(name = "id")]
    ReadId,
    
    /// Reference database mapping: Accession --> ReadKeyword
    #[value(name = "keyword")]
    ReadKeyword,

    /// Reference database mapping: Accession --> ReadAccession
    #[value(name = "accession")]
    ReadAccession,
}