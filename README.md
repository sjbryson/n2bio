<div align="center">
    <img width="480" alt="logo-n2bio 2" src="./assets/n2bio-logo.png" />
</div>

## n2bio - a rust workspace and library for building bioinformatics cli tools

*I created this repo as part of my rust learning journey - building cli tools that I use in my own research & using LLM's along the way.*

---

### Setup & Installation

1. **Install Cargo:** Ensure you have Rust and Cargo installed on your system. If you don't have it, follow the official installation instructions on [rustup.rs](https://rustup.rs/).

2. **Clone & Build:** Clone the repository, navigate into the directory, and build the binary (located in n2bio/target/release/):
   
```
git clone https://github.com/sjbryson/n2bio.git
cd n2bio
cargo build --release
```
3. **Other Useful Commands:**
  - Run Tests: To execute the project's test suite, run: ```cargo test```
  - Local Installation: To install the binary globally to your machine, run: ```cargo install --path . ```
  - For more details on managing Rust packages, visit the official [Cargo Book](https://doc.rust-lang.org/cargo/).

---

### Library:

[**n2bio/**](./n2bio) - Modules I'm developing to work with standard file formats, IO, and common bioinformatics data.
  - sam.rs      - Read and work with SAM formatted alignment records.
  - bam.rs      - Read and work with BAM formatted alignment records.
  - fastq.rs    - Read and write fastq files.
  - fasta.rs    - Read and write fasta files.
  - sequence.rs - Traits to work with DNA sequence data.
  - kmer.rs     - Traits to work with kmers.
  - hist.rs     - Structs and functions to work with distributions and associated stats.
  - metadata.rs - Structs and functions to work with common metadata files (tsv, json, jsonl).
  - readers.rs  - Boilerplate code for reading files and stdin.
  - writers.rs  - Boilerplate code for writing data to files and stdout.

---
### Examples:

[**fastfilter/**](./examples/fastfilter) - Tool to parse SAM formatted stdout from aligners like minimap2, bowtie2, bwa, etc. and write paired reads that pass filter to {prefix}_r1.fq.gz and {prefix}_r2.fq.gz. For use in a pipeline for host read filtering, eliminating some of the common time consuming write-sort-read-filter steps. Unmapped pairs are retained by default. Optional independent alignment quality metrics can also be applied.
- *Example usage for reading and filtering sam records from stdin and writing paired fastq records.*
- *fastfilter is now a subcommand in the [peat cli tool](./peat) as a subcommand* ```peat filter```

[**fastcov/**](./examples/fastcov) - Another tool to parse SAM formatted stdout from aligners like minimap2, bowtie2, bwa, etc. Use in metagenomics pipeline for target identification. Parses SAM records in stdout from aligner, calculates target coverage (per base) and stats. SAM records are passed through to stdout and can be used as input for samtools or written to file. Run and target level stats are writen to .json formatted txt file. All paired primary and secondary alignments that score above at least one set minimum thresholds are writtten to primary and secondary coverage arrays. Mismatch counts are also stored in a mismatch array.
- *Example usage for reading and filtering sam records from stdin.*
- *fastcov is now a subcommand in the [peat cli tool](./peat) as a subcommand* ```peat coverage```

[**pfqbz2gz/**](./examples/pfqbz2gz) - Tool to convert paired fastq records in bz2 format to gz format.
- *Example usage of paired fastq readers and writers.*

---

## Tools Under Development:

<div align="center">
  <p align="center">
    <a href="./peat">
      <img width="180" alt="peat-logo" src="./assets/peat-logo.png" />
    </a>
  </p>
  <p align="center">
    <a href="./peat"><u>P</u>aired-<u>E</u>nd <u>A</u>lignment <u>T</u>ools</a>
  </p>
</div>

#### There are several subcommands for working with paired-end alignment records:

- **peat filter** - Parse SAM records from stdin and filter to create a filtered paired-end fastq.gz library (r1.fq.gz & r2.fq.gz).
- **peat coverage** - Parse SAM records from stdin and calculate coverage for each reference in the sam/bam header
- **peat bam-rep** - Read a name sorted bam file and generate an interactive report
- **peat bin-reads** - Parse SAM records from stdin or BAM and bin read pairs for each target
  
---

**pfqsim/** - Suite of tools to generate synthetic sequencing libraries and test alignment based classification performance.

```
Usage: pfqsim <COMMAND>

Commands:
  model     Build insert size and Q-score distributions from a BAM file
  generate  Generate a simulated paired-read library from a reference FASTA
  compose   Compose a final metagenomic library based on an abundance config
  analyze   Analyze alignments from stdin sam or a bam file
  help      Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version
```

See the [pfqsim README](./pfqsim/README.md) for more information and examples.

  ---