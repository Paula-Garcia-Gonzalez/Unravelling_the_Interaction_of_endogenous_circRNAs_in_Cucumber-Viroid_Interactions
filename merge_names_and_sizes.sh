#!/bin/bash

# Step 1: Extract gene names from the CSV file without the header
tail -n +2 aligned_names.csv > differentiated_names_no_header.txt

# Step 2: Extract genes and sizes from the TAB file without the header
tail -n +2 circular_sizes.tab > circular_sizes_no_header.txt

# Step 3: Sort both files to use `join`
sort differentiated_names_no_header.txt > aligned_names_sorted.txt
sort circular_sizes_no_header.txt > circular_sizes_sorted.txt

# Step 4: Join the files based on common genes and ensure tab-separated output
join -1 1 -2 1 aligned_names_sorted.txt circular_sizes_sorted.txt | awk '{print $1 "\t" $2}' > filtered_with_size.txt

# Step 5: Add header to the output file
echo -e "ID\tLength" | cat - filtered_with_size.txt > filtered_with_size.tab

# Clean up temporary files if needed
rm differentiated_names_no_header.txt circular_sizes_no_header.txt aligned_names_sorted.txt circular_sizes_sorted.txt filtered_with_size.txt

