#!/usr/bin/env bash
# =============================================================================
# Walks-based co-amplification on the hs1 grid, from JaBbA walks to cons80 stats
#
#   SAMPLESHEET_HG19=hg19.csv SAMPLESHEET_HS1=hs1.csv ./run_all.sh          all genes
#   SAMPLESHEET_HG19=hg19.csv SAMPLESHEET_HS1=hs1.csv ./run_all.sh CCND1    one gene
#
# Samplesheet format: see NOTES.md.
#
# Without samplesheets, step 1 is skipped and ARM_DIR must already hold the
# per-build matrices (<ARM_DIR>/<GENE>/{hg19,hs1}/).
#
# Outputs go to ./output/<GENE>/ (override with OUT_DIR), logs to ./logs/.
# =============================================================================
set -euo pipefail
## make every user-supplied path absolute before changing directory
for v in SAMPLESHEET_HG19 SAMPLESHEET_HS1; do
  if [ -n "${!v:-}" ]; then
    [ -f "${!v}" ] || { echo "$v: no such file: ${!v}"; exit 1; }
    export "$v=$(cd "$(dirname "${!v}")" && pwd)/$(basename "${!v}")"
  fi
done
for v in ARM_DIR OUT_DIR; do
  if [ -n "${!v:-}" ]; then mkdir -p "${!v}"; export "$v=$(cd "${!v}" && pwd)"; fi
done
cd "$(dirname "$0")"
export OUT_DIR="${OUT_DIR:-$(pwd)/output}"
export ARM_DIR="${ARM_DIR:-$OUT_DIR/walk_matrices}"
export REF_DIR="${REF_DIR:-$(cd ../../../data/refGenome && pwd)}"
mkdir -p logs
GENES=("$@"); if [ ${#GENES[@]} -eq 0 ]; then GENES=(CCND1 NSD3 EGFR); fi
run () {
  local log="$1"; shift; echo "  -> $log"
  if ! Rscript "$@" > "logs/${log}.log" 2>&1; then
    echo "     FAILED -- see logs/${log}.log"; tail -20 "logs/${log}.log"; exit 1; fi
}
for GENE in "${GENES[@]}"; do
  export GENE
  echo "== $GENE"
  if [ -n "${SAMPLESHEET_HG19:-}" ]; then
    SAMPLESHEET="$SAMPLESHEET_HG19" BUILD=hg19 run "${GENE}_walks_hg19" 01_walks_to_CN_matrix.R
    SAMPLESHEET="${SAMPLESHEET_HS1:?set SAMPLESHEET_HS1}" BUILD=hs1 run "${GENE}_walks_hs1" 01_walks_to_CN_matrix.R
  fi
  run "${GENE}_combine"    02_combine_hg19_hs1.R
  run "${GENE}_plot"       03_plot_coamps.R
  run "${GENE}_stats_chr"  04_BRCAvsNSCLC_FisherStats_cons80.R
  run "${GENE}_stats_lvr"  05_LeftVsRight_FisherStats_cons80.R
done
echo "Done."
