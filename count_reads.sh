#!/bin/bash

#SBATCH --job-name=process_sam_zeros  # Job name to display with squeue
#SBATCH --output=process_sam_zeros_%j.out  # Output file
#SBATCH --ntasks=1  # Number of tasks
#SBATCH --cpus-per-task=2  # Number of CPUs
#SBATCH -t 1-00:00:00  # Runtime in days-hours:minutes:seconds
#SBATCH --qos short  # The QoS to submit the job
#SBATCH --mem-per-cpu=5G  # Memory per CPU in GB

# Directory containing the SAM files
directory="/home/ggomez/tfm_paula/results/quantification"

# Create a directory for temporary files
mkdir -p count_tmp

# Process each SAM file in the directory
for file in "$directory"/*.sam; do
    if [ ! -f "$file" ]; then
        echo "No SAM files found in the directory."
        exit 1
    fi
    name=$(basename "$file" .sam)  # Get the base name without extension
    echo "Processing $file..."

    # Remove duplicates based on columns 1 and 3
    awk '!seen[$1, $3]++' "$file" > "${file}_no_duplicates"
    
    # Extract and count reads, save to a temporary file
    awk '$3 != "*" {print $3}' "${file}_no_duplicates" | sort | uniq -c | awk '{print $2"\t"$1}' > "count_tmp/${name}_counts.txt"
    
    # Sort the temporary file
    sort -k1 "count_tmp/${name}_counts.txt" > "count_tmp/${name}_counts_sorted.txt"
    
    # Delete intermediate no-duplicates file
    rm "${file}_no_duplicates"    
done

# Combine the generated files
echo "Combining counts..."
sorted_files=(count_tmp/*_counts_sorted.txt)
combined_file="/home/ggomez/tfm_paula/results/quantification/final_counts_zeros.txt"

# Use the first file as the base
cp "${sorted_files[0]}" "$combined_file"

columns="0"

# Iterate over the remaining files and merge using `join`, filling missing values with zeros
for ((i=1; i<${#sorted_files[@]}; i++)); do
    current_file="${sorted_files[i]}"
    mv "$combined_file" "${combined_file}.tmp"

    # Add a new column to the output string
    new_column_id=$((i + 1))
    columns+=",1.$new_column_id"
    
    # Add the second column from the second file to the output
    columns_output="$columns,2.2"
    
    echo $columns_output
    
    join -t$'\t' -a1 -a2 -e "0" -o "$columns_output" "${combined_file}.tmp" "$current_file" > "$combined_file"
done

# Clean up temporary files
rm -rf count_tmp "${combined_file}.tmp"

# Add headers to the combined file
echo "Adding headers..."
header="Ref"
for file in "$directory"/*.sam; do
    name=$(basename "$file" .sam)
    header+="\t$name"
done

sed -i "1i $header" "$combined_file"

# Display the result
echo "Combined file ready: $combined_file"


