#!/bin/bash

#SBATCH --job-name=Alignment_Percentage
#SBATCH --output=FileJob_%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH -t 1-00:00:00
#SBATCH --qos short
#SBATCH --mem-per-cpu=5G

input_dir="/home/ggomez/tfm_paula/results/bowtie/sam_output"
output_file="$input_dir/alignment_percentages.tsv"
log_dir="$input_dir/log_alignment"

# Create output file and log directory
mkdir -p "$log_dir"
echo -e "File\tAlignment Percentage" > "$output_file"

# Load required module
module load samtools

# Move to input directory
cd "$input_dir" || exit 1

# Loop over each .sam file and process one by one
for sam in *.sam; do
    base=$(basename "$sam" .sam)
    bam="${base}.bam"
    log_file="${log_dir}/${base}_flagstat.log"

    echo "Processing $sam..."

    samtools view -bS "$sam" > "$bam"
    percentage=$(samtools flagstat "$bam" | tee "$log_file" | awk '/properly paired/ {print $6}')
    echo -e "${base}\t${percentage}" >> "$output_file"
done

echo "Sequential processing complete."
echo "Results saved to: $output_file"
echo "Logs saved in: $log_dir"

