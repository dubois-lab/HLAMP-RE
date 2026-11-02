#!/bin/bash

input_file="samplesheet_tumor_blank.csv"
T2T_ref_dir="../../data/refGenome"

mkdir outputs

while IFS=, read -r col1 col2 col3 col4
do
# Skip the first row (header)
    if $header; then
        header=false
        continue
    fi

    echo "$col1"
cat > jobScripts/job-cnvkit-tumor_${col1}-5kb.sh <<EOF
#!/bin/bash
#SBATCH --job-name=cnvkit-tumor_${col1}-5kb
#SBATCH --time=2-00:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/cnvkit-tumor_${col1}-5kb.out
#SBATCH --error=logs/cnvkit-tumor_${col1}-5kb.err

# compute per bin coverage
cnvkit.py coverage ${col2} ${T2T_ref_dir}/T2T-CHM13v1.1_genomic_UCSC_5kb_bins.bed -o inputs/ -o inputs/${col1}_tumor_5kbBin_UCSC.cnn

# normalize and fix
cnvkit.py fix inputs/${col1}_tumor_5kbBin_UCSC.cnn inputs/${col1}_tumor_5kbBin_UCSC.cnn inputs/CNV_pooled_normal_reference.cnn -o inputs/${col1}_tumor_5kbBin_UCSC.cnr

# segment and call
cnvkit.py segment inputs/${col1}_tumor_5kbBin_UCSC.cnr -o outputs/${col1}_tumor_v_hmmT_5kbBin_UCSC.cns -p 24 -m hmm-tumor -v ${col3}/${col1}.wf_snp.vcf.gz
cnvkit.py call outputs/${col1}_tumor_v_hmmT_5kbBin_UCSC.cns -o outputs/${col1}_tumor_v_hmmT_5kbBin_UCSC.call.cns -v ${col3}/${col1}.wf_snp.vcf.gz --purity ${col4}

EOF
sbatch jobScripts/job-cnvkit-tumor_${col1}-5kb.sh
done < "$input_file"
