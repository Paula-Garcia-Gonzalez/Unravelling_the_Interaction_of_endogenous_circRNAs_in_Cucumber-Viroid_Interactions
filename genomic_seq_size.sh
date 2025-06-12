#!/bin/bash

# Input file
input_file="genomic_sequences.fasta"

# Output file
output_file="circular_sizes.tab"

# Initialize the output file with a header
echo -e "ID\tLength" > "$output_file"

# Variable to store the current ID
current_id=""

# Process the file line by line
while read -r line; do
    if [[ $line == \>* ]]; then
        # If it's a header (line starting with '>'), extract the ID
        current_id=$(echo "$line" | tr -d '>')
    else
        # If it's a sequence, calculate its length and add it to the table
        sequence_length=${#line}
        echo -e "${current_id}\t${sequence_length}" >> "$output_file"
    fi
done < "$input_file"

echo "Base count has been saved to '$output_file'."

