# LinkPrep Hi-C contact maps

These scripts process the LinkPrep Hi-C libraries, whose `.hic` maps are read by the Fig. S4 to S7 scripts as `hic/linkprep/<sample>_1.hic`. All other Hi-C libraries go through nf-core/hic (`../01_run_nfCore_hic.sh`).

| Script | Output in `OUT_DIR` |
|---|---|
| `01_align_and_make_pairs.sh` | `<sample>mapped.pairs`, `<sample>stats.txt`, `<sample>mapped.PT.bam` and its index |
| `02_pairs_to_hic.sh` | `<sample>_1.hic` |

Fill in `samplesheet_linkprep_blank.csv` with one row per sample, then submit from this folder:

```bash
tail -n +2 samplesheet_linkprep_blank.csv | while IFS=, read -r sample fastq_1 fastq_2; do
  sbatch 01_align_and_make_pairs.sh "$sample" "$fastq_1" "$fastq_2"
done
```

When those jobs have finished:

```bash
tail -n +2 samplesheet_linkprep_blank.csv | while IFS=, read -r sample fastq_1 fastq_2; do
  sbatch 02_pairs_to_hic.sh "$sample"
done
```

Without Slurm, run the scripts with `bash` instead of `sbatch`. Set `SLURM_CPUS_PER_TASK` to the number of threads to use.

The scripts need `bwa`, `pairtools` and `samtools` on the `PATH`, Java, and `juicer_tools_1.22.01.jar`. They also need the hs1 reference indexed with `bwa index`, and a tab-separated file of its chromosome names and lengths.

# BWA (https://github.com/lh3/bwa)

Li, H. *Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM*. arXiv, 2013. 1303.3997.

# pairtools (https://github.com/open2c/pairtools)

Open2C et al. *Pairtools: From sequencing data to chromosome contacts*. PLoS Comput Biol, 2024. 20(5): e1012164.

# SAMtools (https://github.com/samtools/samtools)

Danecek, P. et al. *Twelve years of SAMtools and BCFtools*. GigaScience, 2021. 10(2): giab008.

# Juicer tools (https://github.com/aidenlab/juicer)

Durand, N.C. et al. *Juicer provides a one-click system for analyzing loop-resolution Hi-C experiments*. Cell Systems, 2016. 3(1): p. 95-98.
