//! n2bio/pfqsim/src/compose.rs
//! 

use std::io;
use std::fs;
use crate::cli::{ComposeArgs, GenerateArgs};
use crate::config::{Config, Manifest};
use crate::generate;

pub(crate) fn run(args: ComposeArgs) -> io::Result<()> {
    // 1. Load the input configurations
    let config: Config = Config::from_tsv(&args.config)?;

    // 2. Compute the community read distributions
    let manifest: Manifest = Manifest::from_config(&config, args.total_reads, args.abundance_mode);

    // 3. Save the runtime manifest
    let tracking_path: String = format!("{}.manifest.tsv", args.prefix);
    manifest.save_tsv(&tracking_path)?;
    println!("Execution plan mapped and written to: {}", tracking_path);

    // 4. Clean up any pre-existing global target files
    let global_r1: String = format!("{}.r1.fq.gz", args.prefix);
    let global_r2: String = format!("{}.r2.fq.gz", args.prefix);
    let _ = fs::remove_file(&global_r1);
    let _ = fs::remove_file(&global_r2);

    // 5. Execute the manifest actions sequentially, appending sim reads to the global files
    for row in &manifest.rows {
        if row.calculated_reads == 0 { 
            continue; 
        }

        println!("Simulating {} reads for genome: {}", row.calculated_reads, row.id);

        // Map configuration settings into the generate command
        let gen_args: GenerateArgs = GenerateArgs {
            prefix: row.id.clone(), 
            fasta: row.fasta.clone(),
            model: row.model.clone(),
            num_reads: row.calculated_reads,
            read_length: row.read_length,
            sub_rate: row.sub_rate,
            indel_rate: row.indel_rate,
            threads: args.threads,
            circular: row.circular,
            append_mode: true,
            append_path: Some(args.prefix.clone()), 
        };

        generate::run(gen_args)?;
    }

    println!("Simulation complete! {} paired reads saved to:", args.total_reads);
    println!("  R1 -> {}", global_r1);
    println!("  R2 -> {}", global_r2);

    Ok(())
}