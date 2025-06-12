#!/bin/bash

#SBATCH --job-name=FilterBowtieReads      #Job name to show with squeue
#SBATCH --output=filtered_output_%j.out	#Output file
#SBATCH --ntasks=1	# Number of task
#SBATCH --cpus-per-task=2	# Number of cpus
#SBATCH -t 01:00:00	# Runtime in minutes.
#SBATCH --qos short	# The QoS to submit the job.
#SBATCH --mem-per-cpu=5G	# Memory per cpu in G (see also --mem-per-cpu)

# Input and output files
input_dir=/home/ggomez/tfm_paula/results/bowtie_filtrado/aligned
filtered_files_dir=/home/ggomez/tfm_paula/results/bowtie_filtrado/filtered_pairs

# Create temporary directory
mkdir -p tmp
mkdir -p $filtered_files_dir

# Process file pairs
for file1 in "$input_dir"/*_1_paired.fq; do

  # Generate the corresponding paired file name
  file2="${file1/_1_paired.fq/_2_paired.fq}"
  
  file_name1=$(basename "$file1" .fq)
  file_name2=$(basename "$file2" .fq)
  file_name=$(basename "$file1" _1_paired.fq)
  
  # Step 1: Truncate headers using the integrated script
  echo "Truncating FASTQ headers..."
  awk 'NR % 4 == 1 {sub(/ .*/, "", $0)} {print}' "$file1" > tmp/$file_name1"_truncated.fq"
  awk 'NR % 4 == 1 {sub(/ .*/, "", $0)} {print}' "$file2" > tmp/$file_name2"_truncated.fq"
   
  # Step 2: Filter sequences by length using seqkit
  echo "Filtering headers between 1 and 100 nucleotides..."
  seqkit seq -g -m 1 -M 100 -n tmp/$file_name1"_truncated.fq" > tmp/$file_name1"_headers.txt"
  seqkit seq -g -m 1 -M 100 -n tmp/$file_name2"_truncated.fq" > tmp/$file_name2"_headers.txt"
  
  # Concatenate header files and remove duplicates
  cat tmp/$file_name1"_headers.txt" tmp/$file_name2"_headers.txt" | sort | uniq > tmp/$file_name"_concat.txt"
  
  # Add ^ and $ to header names to generate a pattern
  awk '{print "^" $0 "$"}' tmp/$file_name"_concat.txt" > tmp/$file_name"_patterns.txt"
  
  # Filter original files removing invalid sequences
  seqkit grep -rvf tmp/$file_name"_patterns.txt" $file1 > $filtered_files_dir/$file_name1"_filt.fq"
  seqkit grep -rvf tmp/$file_name"_patterns.txt" $file2 > $filtered_files_dir/$file_name2"_filt.fq"
  
  # Remove temporary directory
  #rm -r tmp
  
done

