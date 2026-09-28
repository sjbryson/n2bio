


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