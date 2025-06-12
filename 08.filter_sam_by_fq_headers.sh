#!/bin/bash

#SBATCH --job-name=sam_headers #Job name to show with squeue
#SBATCH --output=sam_headers_%j.out	#Output file
#SBATCH --ntasks=1	# Number of task
#SBATCH --cpus-per-task=2	# Number of cpus
#SBATCH -t 1-00:00:00	# Runtime in minutes.
#SBATCH --qos short	# The QoS to submit the job.
#SBATCH --mem-per-cpu=5G	# Memory per cpu in G (see also --mem-per-cpu)

# Directories
filtered_fq_dir="/home/ggomez/tfm_paula/results/bowtie_filtrado/filtered_pairs"
sam_files_dir="/home/ggomez/tfm_paula/results/bowtie/sam_output"
output_dir="/home/ggomez/tfm_paula/results/quantification"
tmp_dir="${output_dir}/tmp"

# Create output directories if they don't exist
mkdir -p "$output_dir"
mkdir -p "$tmp_dir"

# Process filtered .fq files
for fq_file in "${filtered_fq_dir}"/*1_paired_filt.fq; do
    if [[ -f "$fq_file" ]]; then
        # Get the base name of the file
        base_name=$(basename "$fq_file" _1_paired_filt.fq)

        # Generate file with truncated headers
        headers_file="${tmp_dir}/${base_name}_headers.txt"
        awk 'NR % 4 == 1 {sub(/^@/, "", $0); sub(/ .*/, "", $0)} NR % 4 == 1 {print}' "$fq_file" > "$headers_file"
        echo "Headers extracted for ${fq_file} and saved to ${headers_file}"

        # Filter the corresponding SAM file
        sam_file="${sam_files_dir}/${base_name}.sam"
        if [[ -f "$sam_file" ]]; then
            filtered_sam_file="${output_dir}/${base_name}_filtered.sam"
            grep -F -f "$headers_file" "$sam_file" > "$filtered_sam_file"

            echo "Filtered SAM file saved to ${filtered_sam_file}"
        else
            echo "SAM file not found for ${base_name}, skipping..."
        fi
    else
        echo "No file matching the pattern found in ${filtered_fq_dir}"
    fi
done

echo "Processing completed."

