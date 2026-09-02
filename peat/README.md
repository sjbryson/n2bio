

<p style="text-align:center"><img width="200" alt="peat-logo" src="./assets/peat-logo.png" /></p>

**<p style="text-align: center"><u>P</u>aired-<u>E</u>nd <u>A</u>lignment <u>T</u>ools</p>**

#### There are several subcommands for peat:

```
Usage: peat <COMMAND>

Commands:
  filter     Parse SAM records from stdin and filter to create a filtered paired-end fastq.gz library (r1.fq.gz & r2.fq.gz)
  coverage   Parse SAM records from stdin and calculate coverage for each reference in the sam/bam header
  bam-rep    Read a name sorted bam file and generate an interactive report
  bin-reads  ToDo: Parse SAM records from stdin or BAM and bin read pairs for each target
  help       Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version
  ```

#### peat filter:
```
Usage: peat filter [OPTIONS] --prefix <PREFIX> --report <REPORT> --filter_mode <FILTER_MODE>

Options:
  -t, --threads <THREADS>
          Number of worker threads for parsing and pairing
          
          [default: 4]

      --shards <SHARDS>
          Number of shards for the mate pairing hash map (recommend 4-8x threads, default = 32)
          
          [default: 32]

  -p, --prefix <PREFIX>
          Prefix for output files (e.g. 'out' -> out.r1.fq.gz, out.r2.fq.gz)

  -r, --report <REPORT>
          Name of the run/sample for the JSON report -> creates {report}.json

      --filter_mode <FILTER_MODE>
          How abundance values should be mathematically interpreted

          Possible values:
          - lowpass:  Low Pass: process SAM records that are below defined thresholds
          - highpass: High Pass: process SAM records that are above defined thresholds

      --align_score <ALIGN_SCORE>
          Optional: Alignment Score - sam.get_int_tag("AS")

      --align_length <ALIGN_LENGTH>
          Optional: Alignment Lenth - sam.calculate_alignment_length()

      --base_score <BASE_SCORE>
          Optional: Per base alignment score (AS/AL = avg. align_score per covered base) - sam.calculate_as_al()

      --align_prop <ALIGN_PROP>
          Optional: Alignment Proportion - sam.calculate_alignment_proportion()

      --align_ident <ALIGN_IDENT>
          Optional: Alignment Percent Identity - sam.calculate_alignment_accuracy()

      --mapq <MAPQ>
          Optional: Max MAPQ score - sam.mapq()

  -h, --help
          Print help (see a summary with '-h')
```

#### peat coverage:

```
Usage: peat coverage [OPTIONS] --report <REPORT>

Options:
  -t, --threads <THREADS>            Number of worker threads for parsing [default: 4]
  -r, --report <REPORT>              Name of the run/sample for the JSON report -> creates {report}.json
      --db <DB>                      Optional path to an SQLite taxonomy database (see vref2db)
      --align_score <ALIGN_SCORE>    Optional: Alignment Score - sam.get_int_tag("AS")
      --align_length <ALIGN_LENGTH>  Optional: Alignment Lenth - sam.calculate_alignment_length()
      --base_score <BASE_SCORE>      Optional: Per base alignment score (AS/AL = avg. align_score per covered base) - sam.calculate_as_al()
      --align_prop <ALIGN_PROP>      Optional: Alignment Proportion - sam.calculate_alignment_proportion()
      --align_ident <ALIGN_IDENT>    Optional: Alignment Percent Identity - sam.calculate_alignment_accuracy()
      --mapq <MAPQ>                  Optional: Max MAPQ score - sam.mapq()
  -h, --help                         Print help
```

#### peat bam-rep:

```
Usage: peat bam-rep [OPTIONS] --bam <BAM> --report <REPORT>

Options:
  -b, --bam <BAM>            Input name-sorted BAM file
  -r, --report <REPORT>      Report file prefix - creates {report}.json and optional {report}.html
      --html                 Generate html plots
  -q, --min-mapq <MIN_MAPQ>  Minimum MAPQ score for insert size calculation [default: 40]
  -i, --max-ins <MAX_INS>    Max insert size to use for summary stats calculation [default: 1000]
  -l, --max-len <MAX_LEN>    Max read length to use [default: 150]
  -h, --help                 Print help
  ```

#### peat bin-reads: Under development

```
Usage: peat bin-reads --bam <BAM> --reference-map <REFERENCE_MAP>

Options:
  -b, --bam <BAM>                      Path to an input BAM file to evaluate
  -r, --reference-map <REFERENCE_MAP>  Path to a TSV mapping file: referenc_ id --> bin_id
  -h, --help                           Print help
```

## Examples - 


## Roadmap - 