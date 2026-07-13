---

# pfqsim

A tool to generate synthetic sequencing libraries and analyze read classification

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
    - abundance: [f64]  Relative abundance of total reads to generate for each row's genome. Used to calculate the number of reads to generate for each row - calculation based on proprtion of total reads (default abundance mode "reads") or proportion of genome copies (abundance mode "copies")
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

## Example use:

*Test alignment of simulated reads from 2 reference human genomes and a negative control (mouse genome) to the T2T-CHM13 reference. Example input and output files are located in pfqsim/example_files/*

1. Generate a name sorted bam from a sequencing run that you want to model.

2. Run pfqsim model:

```pfqsim model -b name_sorted.bam -m name_sorted_bam_model```

3. Run pfqsim compose:

```pfqsim compose -c HRFT.config.tsv -m name_sorted_bam_model.json -p HRFT -n 10000000 -t 12 -l 150 --vary-lengths --abundance-mode reads```

4. Align reads to reference genome (e.g. uses minimap2 and a precomputed reference .mmi file).

```minimap2 -ax sr --eqx --secondary=no -t 12 T2T_CHM13.mmi HRFT.r1.fq.gz HRFT.r2.fq.gz | samtools sort -n -o HRFTv13.nsorted.bam -```

5. Run pfqsim analyze:

```pfqsim analyze -c HRFT.manifest.tsv -b HRFTv13.nsorted.bam -r CHM13_mapping.tsv -o HRFTv13_analyze -m keyword```

6. Optionally run bamrep to examine the bam file generated by the alignment.

```bamrep -b HRFTv13.nsorted.bam -r HRFTv13_bamrep --html```
