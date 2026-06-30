//! n2bio/pfqsim/src/compose.rs
//! 

use std::io;
use std::path::PathBuf;

use crate::cli::{ ComposeArgs, GenerateArgs };
use crate::config::{ Config, Manifest };
use crate::generate;

pub(crate) fn run(args: ComposeArgs) -> io::Result<()> {
    // 1. Load the input configurations
    let config: Config = Config::from_tsv(&args.config)?;

    // 2. Compute the community read distributions
    let manifest: Manifest = Manifest::from_config(&config, args.total_reads, args.abundance_mode);

    // 3. Save the runtime manifest
    let tracking_path: String = format!("{}_manifest.tsv", args.prefix.display());
    manifest.save_tsv(&tracking_path)?;
    println!("Execution plan mapped and written to: {}", tracking_path);

    // 4. Execute the manifest actions sequentially
    for row in &manifest.rows {
        if row.calculated_reads == 0 { continue; }

        let (r1_path, r2_path) = match (&row.r1_fq, &row.r2_fq) {
            (Some(r1), Some(r2)) if r1.exists() && r2.exists() => {
                (r1.clone(), r2.clone())
            }
            _ => {
                let out_prefix: String = format!("{}_sim", row.id);
                let gen_args: GenerateArgs = GenerateArgs {
                    fasta: row.fasta.clone(),
                    model: row.model.clone(),
                    num_reads: row.calculated_reads,
                    read_length: row.read_length,                 // read_length
                    sub_rate: row.sub_rate,
                    indel_rate: row.indel_rate,
                    threads: args.threads,
                    prefix: out_prefix.clone(),
                    circular: row.circular,
                };
                generate::run(gen_args)?;
                
                (PathBuf::from(format!("{}.r1.fq.gz", out_prefix)), 
                 PathBuf::from(format!("{}.r2.fq.gz", out_prefix)))
            }
        };
        
        // TODO Step 5: Queue these paths into your upcoming master mixer/shuffler loop!
    }

    Ok(())
}