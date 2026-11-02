#!/bin/bash

input_file="samplesheet_blank.csv"

while IFS=, read -r col1
do
    if $header; then
        header=false
        continue
    fi
    echo "Processing: $col1"

cat > walks/jobScripts/job-script_${col1}_jabba_walks.sh <<EOF
#!/bin/bash
#SBATCH --job-name=${col1}-jabba_walks
#SBATCH --output=logs/jabba_walks_${col1}.out
#SBATCH --error=logs/jabba_walks_${col1}.err
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=5
#SBATCH --time=1-12:00:00
#SBATCH --mem=40G

Rscript 01_JaBbA_walks.R $col1

EOF
sbatch walks/jobScripts/job-script_${col1}_jabba_walks.sh
done < "$input_file"
