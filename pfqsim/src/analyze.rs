//! n2bio/pfqsim/src/compose.rs
//! 
use std::io;

use crate::cli::AnalyzeArgs;


pub fn run(args: AnalyzeArgs) -> io::Result<()> {
    println!("Analyzing alignments from bam file: {:?}", args.bam);
    
    Ok(())
}