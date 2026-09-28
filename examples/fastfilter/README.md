


**fastfilter/** - Tool to parse SAM formatted stdout from aligners like minimap2, bowtie2, bwa, etc. and write paired reads that pass filter to {prefix}_r1.fq.gz and {prefix}_r2.fq.gz. For use in a pipeline for host read filtering, eliminating some of the common time consuming write-sort-read-filter steps. Unmapped pairs are retained by default. Optional independent alignment quality metrics can also be applied.

**Pipeline example:**

```
minimap2 -ax sr --eqx --secondary=no {map_threads} {input_mmi} {r1} {r2} | \
fastfilter {filter_threads} {max_ap} {max_pi} {max_as} {max_al} {max_sl} {max_mq} {fq_prefix}


Usage: fastfilter [OPTIONS] --fq-prefix <FQ_PREFIX>

Options:
  -t, --threads <THREADS>      Number of worker threads for parsing and pairing [default: 4]
      --shards <SHARDS>        Number of shards for the ShardedMateMap (recommend 4-8x threads) [default: 32]
  -p, --fq-prefix <FQ_PREFIX>  Prefix for output files (e.g. 'out' -> out_r1.fq.gz, out_r2.fq.gz)
      --max-ap <MAX_AP>        Optional: Max Alignment Proportion
      --max-ai <MAX_AI>        Optional: Max Alignment Identity
      --max-as <MAX_AS>        Optional: Max Alignment Score
      --max-al <MAX_AL>        Optional: Max Alignment Lenth
      --max-bs <MAX_BS>        Optional: Max per base alignment score
      --max-mq <MAX_MQ>        Optional: Max MAPQ score
  -h, --help                   Print help
  -V, --version                Print version
```