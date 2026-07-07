//! n2bio/pfqsim/src/compose.rs
//! 

use std::io;
use std::fs;
use std::time::Instant;

use crate::cli::{ComposeArgs, GenerateArgs};
use crate::config::{ComposeConfig, Manifest};
use crate::generate;

pub(crate) fn run(args: ComposeArgs) -> io::Result<()> {
    let start_time: Instant = Instant::now();
    // 1. Load the input configuration, calculate genome lengths, and validate circularity
    let mut config: ComposeConfig = ComposeConfig::from_tsv(&args.config)?;
    println!("Validating reference FASTA metrics and parsing genome lengths...");
    config.validate_and_compute_lengths()?;
    
    // 2. Compute the community read distributions
    let manifest: Manifest = Manifest::from_config(&config, args.total_reads, args.abundance_mode);

    // 3. Save the runtime manifest
    let manifest_path: String = format!("{}.manifest.tsv", args.prefix);
    manifest.save_tsv(&manifest_path)?;
    println!("Library manifest written to: {}", manifest_path);

    // 4. Overwrite behavior - delete pre-existing r1 and r2 fastq
    //    Could change to raise overwrite warning...
    let global_r1: String = format!("{}.r1.fq.gz", args.prefix);
    let global_r2: String = format!("{}.r2.fq.gz", args.prefix);
    fs::remove_file(&global_r1).ok();
    fs::remove_file(&global_r2).ok();

    // 5. Execute the manifest actions sequentially, appending sim reads to the global files
    for row in &manifest.rows {
        println!("Simulating {} reads for genome: {}", row.calculated_reads, row.id);
        
        if row.calculated_reads == 0 { 
            continue; 
        }
        
        // Map configuration settings into the generate command
        let gen_args: GenerateArgs = GenerateArgs {
            prefix: row.id.clone(), 
            keyword: row.keyword.clone(),
            fasta: row.fasta.clone(),
            model: args.model.clone(),
            num_reads: row.calculated_reads,
            read_length: args.read_length,
            sub_rate: row.sub_rate,
            indel_rate: row.indel_rate,
            threads: args.threads,
            circular: row.circular,
            vary_lengths: args.vary_lengths,
            append_mode: true,
            append_path: Some(args.prefix.clone()), 
        };

        generate::run(gen_args)?;
    }

    println!("Simulation complete! {} paired reads saved to:", args.total_reads);
    println!("  R1 -> {}", global_r1);
    println!("  R2 -> {}", global_r2);
    println!("Runtime: {:.2} seconds", start_time.elapsed().as_secs_f64());

    Ok(())
}