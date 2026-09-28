


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

An optional metadata file (--metadata or -m) can be used to add additional information for each reference sequence in the coverage report. The --metadata_key or -k option tells fastcov which column or field in the metadata file corresponds to the reference sequence identifier - e.g. a column named "accession" could refer to the accessions in the reference database that was aligned to - these should match what you would see in a sam/bam header. All additional fields and values associated with each key will be included in the report.json file.

```
Usage: fastcov [OPTIONS] --run-name <RUN_NAME>

Options:
  -t, --threads <THREADS>    Number of worker threads for parsing and pairing [default: 4]
  -r, --report <REPORT>      Name of the run/sample for the JSON report -> creates {report}.json
  -m, --metadata <METADATA>  Optional path to a metadata file
  -k, --metadata-key <METADATA_KEY>  Optional metadata keyword
      --min-ap <MIN_AP>      Optional: Min Alignment Proportion
      --min-ai <MIN_AI>      Optional: Min Alignment Identity
      --min-as <MIN_AS>      Optional: Min Alignment Score
      --min-al <MIN_AL>      Optional: Min Alignment Lenth
      --min-bs <MIN_BS>      Optional: Min per base alignment score
      --min-mq <MIN_MQ>      Optional: Min MAPQ score
  
  -h, --help                 Print help
  -V, --version              Print version
```