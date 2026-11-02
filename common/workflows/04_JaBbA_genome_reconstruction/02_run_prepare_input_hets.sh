#!/bin/bash

input_file="samplesheet_tumor_blank.csv"

while IFS=, read -r col1 col2 col3
do
# Skip the first row (header)
    if $header; then
        header=false
        continue
    fi
    
cat > inputdata/hets/intermediate_results/job-script_${col1}_pysam_vcftohets.sh <<EOF
#!/bin/bash
#SBATCH --job-name=${col1}-pysam_vcftohets
#SBATCH --output=inputdata/hets/intermediate_results/${col1}-pysam_vcftohets.out
#SBATCH --error=inputdata/hets/intermediate_results/${col1}-pysam_vcftohets.err
#SBATCH --mem=156G
#SBATCH --cpus-per-task=48
#SBATCH --time=1-12:00

python 02_prepare_input_hets.py ${col1} ${col3}

# remove chr prefix
awk 'BEGIN{OFS="\t"} {$1 = gensub(/^chr/, "", 1, $1); print}' inputdata/hets/intermediate_results/${col1}_hets_pysam.tsv > inputdata/hets/${col1}_hets_pysam_noChr.tsv

EOF
sbatch inputdata/hets/intermediate_results/job-script_${col1}_pysam_vcftohets.sh
done < "$input_file"
