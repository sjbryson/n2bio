//! n2bio/pfqsim/src/main.rs
//! 

mod cli;
mod model;
mod generate;
mod compose;
mod analyze;
mod inserts;
mod qualities;
mod simstats;
mod genome;
mod mutate;

use std::io;
use clap::Parser;
use cli::{ Cli, Commands };


fn main() -> io::Result<()> {
    // Parse the command line arguments
    let args = Cli::parse();

    // Match on the subcommand and route to the correct run function
    match args.command {
        Commands::Model(model_args)       => model::run(model_args)?,
        Commands::Generate(gen_args)   => generate::run(gen_args)?,
        Commands::Compose(comp_args)    => compose::run(comp_args)?,
        Commands::Analyze(analyze_args) => analyze::run(analyze_args)?,
    }

    Ok(())
}