//! n2bio/peat/src/main.rs
//! 

#![allow(unused)]

mod cli;
mod filter;
mod coverage;
mod bamrep;
mod binreads;
mod samfilters;
mod bamrep_stats;
mod bamrep_report;

use std::io;
use clap::Parser;
use crate::cli::{ Cli, Commands };


fn main() -> io::Result<()> {
    // Parse the command line arguments
    let args: Cli = Cli::parse();

    // Match on the subcommand and route to the correct run function
    match args.command {
        Commands::Filter(filter_args)       => filter::run(filter_args)?,
        Commands::Coverage(coverage_args) => coverage::run(coverage_args)?,
        Commands::BamRep(bamrep_args)       => bamrep::run(bamrep_args)?,
        Commands::BinReads(binreads_args) => binreads::run(binreads_args)?,
    }

    Ok(())
}