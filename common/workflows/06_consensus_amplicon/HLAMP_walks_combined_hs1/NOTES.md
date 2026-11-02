# Walks-based consensus amplicon, hg19 + hs1 combined

Co-amplification frequencies around CCND1 (chr11), NSD3 (chr8) and EGFR (chr7)
from JaBbA walks, for BRCA and NSCLC (LUSC), on the hs1 (T2T-CHM13 v2.0) 10 kb
grid, and the cons80 Fisher tests built on them. Starts from the JaBbA walk
files (gGnome `gWalk` `.rds`, see `05_JaBbA_genome_walks`).

## Input

Two samplesheets, one row per walk file (templates:
`samplesheet_walks_{hg19,hs1}_blank.csv`):

| column | |
|---|---|
| `sample` | sample ID; becomes the matrix column |
| `cohort` | `BRCA` / `NSCLC` (hg19), `BRCA` / `LUSC` (hs1) |
| `patient` | patient ID; at most one sample per patient is kept (empty = no dedup) |
| `walk_file` | path to the walk `.rds` |
| `priority` | tie-break between walk sources, lower wins (optional) |

* **hg19** — HMF and PCAWG samples, whose JaBbA graphs and walks are on hg19.
  The analysis used 10,103 walk files from 631 samples. Cohort is BRCA for
  breast and NSCLC for lung primaries. `patient` is the HMF patient ID, so
  only one biopsy per HMF patient is kept.
* **hs1** — ONT samples reconstructed on hs1. The analysis used 982 walk files
  from 56 samples; the same samplesheet serves all three genes.

The samplesheets and walks contain patient-level data and are not part of
this repository.

Reference files, from `common/data/refGenome/`: `hg19ToHs1.over.chain`,
`hs1.chrom.sizes`, `cgc_T2T_gr_consAmpWalks.rds`, `T2T_gencode_gt_UCSC.rds`.

## Steps

| # | script | what |
|---|---|---|
| 1 | `01_walks_to_CN_matrix.R` | per build: qualifying walks → collapsed segments → per-sample CN matrix (run with `BUILD=hg19` and `BUILD=hs1`) |
| 2 | `02_combine_hg19_hs1.R` | project the hg19 matrix onto hs1 and merge with the hs1 matrix |
| 3 | `03_plot_coamps.R` | per-cohort profiles (%amp, AmpCount, MedianCN): gene ±5 Mb, ±2.5 Mb, chromosome, genome |
| 4 | `04_BRCAvsNSCLC_FisherStats_cons80.R` | BRCA vs NSCLC Fisher test per bin on the gene chromosome |
| 5 | `05_LeftVsRight_FisherStats_cons80.R` | left-vs-right mirror-pair Fisher test per cohort, gene ±2.5 Mb |

```bash
SAMPLESHEET_HG19=/path/hg19.csv SAMPLESHEET_HS1=/path/hs1.csv ./run_all.sh
```

Requires R with gGnome, gUtils, gTrack, data.table, GenomicRanges, rtracklayer,
BSgenome.Hsapiens.UCSC.hg19, qvalue (matrixStats optional). Step 1 on the hg19
walks reads every walk file once per gene; set `MC_CORES` to parallelise.

## Method

**Qualifying walks and samples (step 1).** A walk qualifies if its segments at
CN ≥ 10 cover ≥ 30 % of the gene's target interval:

| gene | hg19 | hs1 |
|---|---|---|
| CCND1 | 11:69,455,873–69,469,242 (NM_053056 coding sequence) | chr11:69,658,031–69,671,351 |
| NSD3 | 8:38,127,215–38,239,790 (CGC, WHSC1L1) | chr8:38,546,659–38,659,240 |
| EGFR | 7:55,086,714–55,324,313 (CGC) | chr7:55,178,937–55,372,056 |

hs1 intervals are the T2T gene elements from `cgc_T2T_gr_consAmpWalks.rds`.

All qualifying walks of a sample, from all of its walk files, are collapsed
into one segment set. A sample with at least one qualifying walk enters the
matrix. Where samples share a `patient`, the one with the highest gene coverage
is kept. Ties go to the higher CN at the gene of its best walk, then to the
sample ID. Each 10 kb bin of the build's genome gets, per sample, the maximum CN
over the collapsed segments hitting it. A bin no segment of the sample hits
gets CN 2.

**Projection (step 2).** Every hs1 10 kb tile inherits its CN from the hg19
10 kb bin that covers most of it after liftOver (`hg19ToHs1.over.chain`). The
winning bin is chosen against the complete hg19 10 kb grid
(`tileGenome(seqlengths(BSgenome.Hsapiens.UCSC.hg19), 10000)`), and only then
looked up in the hg19 matrix. If the winning bin carries no walk for a sample,
the cell is `NA`. hg19 bins that do not lift, or lift onto more than one
chromosome, are dropped. Rows are the tiles of the hs1 matrix plus the hs1
tiles spanned by the lifted hg19 walk bins.

**Per-bin summaries.** A sample is amplified in a bin at CN ≥ 10.
`AmpSum` = number of amplified samples, `AmpPercent` = `AmpSum / n_samples_total`
(all samples of the cohort), `n_samples_obs` = samples with a non-`NA` cell.

**cons80 tests (steps 4, 5).** Only bins (or mirror pairs) where at least one
side reaches an amplification frequency ≥ 80 % are tested, with a 2×2 Fisher
exact test of amplified vs not amplified. q-values by `qvalue()`, with
Benjamini–Hochberg as fallback.

## Output

```
output/
├── walk_matrices/<GENE>/{hg19,hs1}/        step 1: per-build CN matrices, collapsed segments, sample_status.csv
└── <GENE>/
    ├── {BRCA,NSCLC_LUSC}_<GENE>_CNMatm_combined_hs1.{txt.gz,rds}
    ├── {BRCA,NSCLC_LUSC}_<GENE>_CNMatm_combined_hs1_coamps.pdf
    └── stats/cons80/
        ├── <GENE>_BRCAvsNSCLC_FisherStats_CN10_cons80.txt            all tested bins on the gene chromosome
        ├── <GENE>_BRCAvsNSCLC_25e5window_FisherStats_CN10_cons80.txt  subset in gene ±2.5 Mb
        ├── <GENE>_BRCAvsNSCLC_FisherStats_CN10_cons80_gr.rds
        ├── <GENE>_BRCAvsNSCLC_FisherStats_CN10_cons80_p0.05.pdf
        ├── <GENE>_{BRCA,NSCLC}_LeftVsRight_FisherStats_CN10_cons80.txt
        └── <GENE>_BRCA_NSCLC_LeftVsRight_FisherStats_CN10_cons80_p0.05.pdf
```

`<GENE>_BRCAvsNSCLC_FisherStats_CN10_cons80_gr.rds` and
`<GENE>_{BRCA,NSCLC}_LeftVsRight_FisherStats_CN10_cons80.txt` are the files
read by the Figure 2E / 3A / 3B scripts from `IN_DIR`.
