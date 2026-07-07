<img width="3390" height="1200" alt="logo-n2bio 2" src="https://github.com/user-attachments/assets/2adfdafa-0ba3-4aa8-a891-75bb35e3871f" />

## n2bio - a rust workspace for building bioinformatics cli tools

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

### Modules:

**n2core/** - Modules I'm developing to work with standard file formats and IO.

---

### Cli Tools:

**fastfilter/** - Tool to parse SAM formatted stdout from aligners like minimap2, bowtie2, bwa, etc. and write paired reads that pass filter to {prefix}_r1.fq.gz and {prefix}_r2.fq.gz. For use in a pipeline for host read filtering, eliminating some of the common time consuming write-sort-read-filter steps. Unmapped pairs are retained by default. Optional independent alignment quality metrics can also be applied.

**Pipeline example:**

```
minimap2 -ax sr --eqx --secondary=no {map_threads} {input_mmi} {r1} {r2} | \
fastfilter {filter_threads} {max_ap} {max_pi} {max_as} {max_al} {max_sl} {max_mq} {fq_prefix}


Usage: fastfilter [OPTIONS] --fq-prefix <FQ_PREFIX>

Options:
  -t, --threads <THREADS>      Number of worker threads for parsing and pairing [default: 4]
      --shards <SHARDS>        Number of shards for the ShardedMateMap (recommend 4-8x threads) [default: 64]
  -p, --fq-prefix <FQ_PREFIX>  Prefix for output files (e.g. 'out' -> out_r1.fq.gz, out_r2.fq.gz)
      --max-ap <MAX_AP>        Optional: Max Alignment Proportion
      --max-pi <MAX_PI>        Optional: Max Percent Identity
      --max-as <MAX_AS>        Optional: Max Alignment Score
      --max-al <MAX_AL>        Optional: Max Alignment Lenth
      --max-sl <MAX_SL>        Optional: Max AS/AL score
      --max-mq <MAX_MQ>        Optional: Max MAPQ score
  -h, --help                   Print help
  -V, --version                Print version
```

---

**fastcov/** - Another tool to parse SAM formatted stdout from aligners like minimap2, bowtie2, bwa, etc. Use in metagenomics pipeline for target identification. Parses SAM records in stdout from aligner, calculates target coverage (per base) and stats. SAM records are passed through to stdout and can be used as input for samtools or written to file. Run and target level stats are writen to .json formatted txt file. All paired primary and secondary alignments that score above at least one set minimum thresholds are writtten to primary and secondary coverage arrays. Mismatch counts are also stored in a mismatch array.

**Pipeline example:**

```
minimap2 -ax sr --eqx {map_threads} {input_mmi} {r1} {r2} | \
fastcov {cov_threads} -r {sample} {min_as} | \
samtools sort {sort_threads} - -o {sample}.sorted.bam
```
Or if you don't want to save the sam/bam file - pipe to /dev/null:
```
minimap2 -ax sr --eqx {map_threads} {input_mmi} {r1} {r2} | \
fastcov {cov_threads} -r {sample} {min_as} > /dev/null
```
And if you want to test filtering parameters from an existing sam/bam file:
```
samtools view -h file.bam | fastcov {cov_threads} -r {sample} {min_as} > /dev/null
```
If a viral taxonomy db was created using vref2db (use option --db <path to SQLite db file>) lineage data for the associated Accession will be reported.
```
Usage: fastcov [OPTIONS] --run-name <RUN_NAME>

Options:
  -t, --threads <THREADS>    Number of worker threads for parsing and pairing [default: 4]
  -r, --run-name <RUN_NAME>  Name of the run/sample for the JSON report
      --min-ap <MIN_AP>      Optional: Min Alignment Proportion
      --min-pi <MIN_PI>      Optional: Min Percent Identity
      --min-as <MIN_AS>      Optional: Min Alignment Score
      --min-al <MIN_AL>      Optional: Min Alignment Lenth
      --min-sl <MIN_SL>      Optional: Min AS/AL score
      --min-mq <MIN_MQ>      Optional: Min MAPQ score
      --db <DB>              Optional path to an SQLite taxonomy database (see vref2db)
  -h, --help                 Print help
  -V, --version              Print version
```

---

**vref2db/** - NCBI Virus-Host taxonomy DB builder. Create an SQLite database for a set of reference viral sequences (NCBI Accessions). 

1. Start by generating a list of Accession ID's from a fasta file:

```
grep ">" viral_refs.fna | awk '{print $1}' | sed 's/>//' > viral_refs.accessions.txt
```

2. Use NCBI's **datasets** to download taxonomic information in json format.

```
datasets summary virus genome accession --inputfile viral_refs.accessions.txt > viral_refs.metadata.tsv
```

3. Download the ncbi [taxdump](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/) to get the nodes.dmp file.

4. Run **vref2db**:

```
vref2db --db db_name --json viral_refs.metadata.json --nodes nodes.dmp

Usage: vref2db --db <DB> --json <JSON> --nodes <NODES>

Options:
      --db <DB>        Name of the SQLite database to create (e.g., taxonomy.db)
      --json <JSON>    Path to the JSON report generated by NCBI Datasets
      --nodes <NODES>  Path to the nodes.dmp file
  -h, --help           Print help
  -V, --version        Print version

```
---

**bamrep/** - Generate a report from a name sorted bam file for an alignment of paired end reads. Optional Html and Svg reports in addition to standard json output of summary stats + histograms for insert size (calculated using alignments that pass user defined max insert size and min mapq values), mapq values, alignment scores, alignment lengths, per base alignments cores, alignment proportions, and alignment percent identities. Html output (--html option) uses plotly, has sliders to check thresholds for filtering alignments, and allows savingindividual plots as svg files.

```
Usage: bamrep [OPTIONS] --bam <BAM> --report <REPORT>

Options:
  -b, --bam <BAM>            Input name-sorted BAM file
  -r, --report <REPORT>      Output JSON report file
      --html                 Generate html plots
  -q, --min-mapq <MIN_MAPQ>  Minimum MAPQ score for insert size calculation [default: 40]
  -i, --max-ins <MAX_INS>    Max insert size to use for summary stats calculation [default: 1000]
  -l, --max-len <MAX_LEN>    Max read length to use [default: 150]
  -h, --help                 Print help
  -V, --version              Print version
```

---

**pfqbz2gz/** - Tool to convert paired fastq records in bz2 format to gz format.

```
Usage: pfqbz2gz [OPTIONS] --r1 <R1> --r2 <R2> --output-prefix <OUTPUT_PREFIX>

Options:
  -1, --r1 <R1>                        Path to R1 bz2 file
  -2, --r2 <R2>                        Path to R2 bz2 file
  -o, --output-prefix <OUTPUT_PREFIX>  Output prefix for the new gz files (e.g. 'sample1' becomes 'sample1_R1.fq.gz')
  -t, --threads <THREADS>              Total CPU threads to allocate across the pipeline [default: 4]
  -h, --help                           Print help
  -V, --version                        Print version
```

---

**pfqsim/** - Tool to generate synthetic sequencing libraries.

```
Fast metagenomic read simulator

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

- **pfqsim model** - Build insert size, read length, and Q-score distributions from a name sorted BAM file. Make sure the bam was name sorted e.g. samtools sort -n -o name_sorted.bam -

```
Usage: pfqsim model [OPTIONS] --bam <BAM> --model <MODEL>

Options:
  -b, --bam <BAM>                  Path to the BAM file for modeling insert size, read length, and Qscore distributions
  -m, --model <MODEL>              Name for the JSON model report -> creates {model}.json
  -l, --read-length <READ_LENGTH>  Read length to model [default: 150]
  -q, --mapq <MAPQ>                Min MAPQ score for filtering alignments for insert size distribution [default: 40]
  -i, --max-ins <MAX_INS>          Max insert size to use for insert size distribution [default: 1000]
  -h, --help                       Print help
```
- **pfqsim generate** - Generate a simulated paired-read library from a reference FASTA
```
Usage: pfqsim generate [OPTIONS] --model <MODEL> --fasta <FASTA> --num-reads <NUM_READS> --prefix <PREFIX> --keyword <KEYWORD> --threads <THREADS>

Options:
  -m, --model <MODEL>              Path to the JSON model report -> created by pfqsim model
  -f, --fasta <FASTA>              Path to the genome fasta file to generate reads from
  -c, --circular                   Boolean value: circularize genome before generating reads
  -s, --sub-rate <SUB_RATE>        Float value for random substitution rate to apply to simulated reads (range: 0.0 - 1.0) [default: 0]
  -i, --indel-rate <INDEL_RATE>    Float value for random insertion and deletion rate to apply to simulated reads (range: 0.0 - 1.0) [default: 0]
  -n, --num-reads <NUM_READS>      Integer value for number of paired reads to create (1 = 1 R1.fq.gz + 1 R2.fq.gz)
  -l, --read-length <READ_LENGTH>  Read length to generate [default: 150]
      --vary-lengths               Boolean value: Vary read lengths based on model
  -p, --prefix <PREFIX>            Prefix for output fastq.gz files (e.g. {prefix}.r1.fq.gz) and for read identifiers (e.g. @{prefix}:{keyword}:Accession::Read Num)
  -k, --keyword <KEYWORD>          Additional keyword to add to read identifiers for use in query-target mapping (e.g. @{prefix}:{keyword}:Accession::Read Num)
  -t, --threads <THREADS>          Number of worker threads
  -h, --help                       Print help
  ```
- **pfqsim compose** - Compose a final metagenomic library based on an abundance config. 
```
Usage: pfqsim compose [OPTIONS] --config <CONFIG> --model <MODEL> --prefix <PREFIX> --total-reads <TOTAL_READS> --threads <THREADS>

Options:
  -c, --config <CONFIG>
          Path to a TSV config file

  -m, --model <MODEL>
          Path to the JSON model report -> created by pfqsim model

  -p, --prefix <PREFIX>
          Prefix for the manifest tsv and both simulated reads (R1 & R2) files

  -n, --total-reads <TOTAL_READS>
          Integer value for number of paired reads to create (1 = 1 R1.fq.gz + 1 R2.fq.gz)

  -l, --read-length <READ_LENGTH>
          Read length to generate
          
          [default: 150]

      --vary-lengths
          Boolean value: Vary read lengths based on model

  -t, --threads <THREADS>
          Number of worker threads

      --abundance-mode <ABUNDANCE_MODE>
          How abundance values should be mathematically interpreted

          Possible values:
          - reads:  Calculate abundance as fraction of total reads
          - copies: Calculate abundance as fraction of total genome copies
          
          [default: reads]

  -h, --help
          Print help (see a summary with '-h')
```

The configuration file requires the following fields:

    - id: [String]  This will be incorporated into the sequence identifier for all reads generated from each row's genome ("@{id}:...).
    - keyword: [String]  This is incroporated as the second position in each sequence identifier (e.g. @{id}:{keyword}:...) 
    - abundance: [f64]  Relative abundance of total reads (range: 0-1) to generate for each row's genome. Used to calculate the number of reads to generate for each row - calculation based on proprtion of total reads (default abundance mode "reads") or proportion of genome copies (abundance mode "copies")
    - fasta: [Path] Path to the genome fasta for this row.
    - circular: [bool]  Should the genome in this row be circularized before read generation. Only works with single contig genomes.
    - sub_rate: [f64] Substitution rate (range: 0-1) to apply to simulated reads.
    - indel_rate: [f64] Random insertion or deletion rate (range: 0-1) to apply to simulated reads.
    
A validation step calculates all genome sizes for abundance calculations and forces circular=false for any multi contig genome fasta.
<br>

- **pfqsim analyze** - Analyze alignments from a name sorted BAM file for classification stats, e.g. accuracy, sensitivity, ROC, etc. Generates an html report with interactive threshholding and json report.

```
Usage: pfqsim analyze [OPTIONS] --config <CONFIG> --bam <BAM> --reference-map <REFERENCE_MAP> --mapping-mode <MAPPING_MODE>

Options:
  -c, --config <CONFIG>
          Path to a TSV config file

  -b, --bam <BAM>
          Path to an input BAM file to evaluate

  -r, --reference-map <REFERENCE_MAP>
          Path to a TSV mapping file: reference id --> mapping-mode(id, keyword, or accession)

  -o, --output <OUTPUT>
          Output path prefix for the generated HTML & json evaluation reports
          
          [default: pfqsim_report]

  -m, --mapping-mode <MAPPING_MODE>
          Part of read identifier to use for reference sequence mapping

          Possible values:
          - id:        Reference database mapping: Accession --> ReadId
          - keyword:   Reference database mapping: Accession --> ReadKeyword
          - accession: Reference database mapping: Accession --> ReadAccession

  -h, --help
          Print help (see a summary with '-h')
```

  ---

**xbag/** - Very experimental code for building nucleotide variation graphs from a set of closely related (e.g. strain/species) genomes.