#!/bin/bash

input_file="samplesheet_tumor_blank.csv"

while IFS=, read -r col1 col2 col3
do
# Skip the first row (header)
    if $header; then
        header=false
        continue
    fi
    
    mkdir -p results/${col1}_JaBbA/logs
    touch results/${col1}_JaBbA/logs/logs_${col1}_jabbaRun.txt
    cd results/${col1}_JaBbA

cat > job-script_${col1}_jabba.sh <<EOF
#!/bin/bash
#SBATCH --job-name=${col1}-jabba
#SBATCH --output=logs/${col1}-jabba.out
#SBATCH --error=logs/${col1}-jabba.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --time=1-01:00
#SBATCH --mem=256G

CPXPARAM_Threads=24
echo "Processing sample: ${col1}"
jba inputdata/SV_bedpe/${col1}_sniffles_PoNfiltered.bedpe inputdata/coverage/${col1}_tumor_5kbBin_UCSC_cov.csv \
    -o results/${col1}_JaBbA \
    --purity ${col2} -n ${col1} --cores 24 -v --mem 256 --ism FALSE --ppmethod='ppgrid' --slack 100 \
    --ploidy 1.5,4.2 \
    --tilim 43200 \
    --hets inputdata/hets/${col1}_hets_pysam_noChr.tsv --seg inputdata/segmentation/hmm/${col1}_tumor_v_hmmT_5kbBin_UCSC.rds \
    > logs/logs_${col1}_jabbaRun.txt

EOF
sbatch job-script_${col1}_jabba.sh
done < "$input_file"

