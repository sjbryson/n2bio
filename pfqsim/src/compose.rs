//! n2bio/pfqsim/src/compose.rs
//! 
use std::io;

use crate::cli::ComposeArgs;


pub fn run(args: ComposeArgs) -> io::Result<()> {
    println!("Generating library from {:?}", args.config);
    
    Ok(())
}