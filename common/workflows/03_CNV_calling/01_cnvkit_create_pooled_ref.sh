#!/bin/bash

input_file="samplesheet_normals_blank.csv"
T2T_ref_dir="../../data/refGenome"

mkdir jobScripts logs inputs

while IFS=, read -r col1 col2
do
# Skip the first row (header)
    if $header; then
        header=false
        continue
    fi

    echo "$col1"
cat > jobScripts/job-cnvkit-normals_${col1}_5kbBin.sh <<EOF
#!/bin/bash
#SBATCH --job-name=cnvkit-normals-${col1}_5kb
#SBATCH --time=0-20:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/cnvkit-normals-${col1}_5kb.out
#SBATCH --error=logs/cnvkit-normals-${col1}_5kb.err

# compute per bin coverage
cnvkit.py coverage ${col2} ${T2T_ref_dir}/T2T-CHM13v1.1_genomic_UCSC_5kb_bins.bed -o inputs/${col1}_normal_5kbBin_UCSC.cnn

# build pooled reference
cnvkit.py reference inputs/*_normal_5kbBin_UCSC.cnn -f ${T2T_ref_dir}/UCSC_T2T-CHM13v1.1_genomic.fna -o inputs/CNV_pooled_normal_reference.cnn

EOF
sbatch jobScripts/job-cnvkit-normals_${col1}_5kbBin.sh
done < "$input_file"
