#!/bin/bash

#SBATCH --job-name=hisat2_align    #Job name to show with squeue
#SBATCH --output=hisat2_align_%j.out	#Output file
#SBATCH --ntasks=16	# Number of task
#SBATCH --cpus-per-task=4	# Number of cpus
#SBATCH -t 1-00:00:00	# Runtime in minutes.
#SBATCH --qos short	# The QoS to submit the job.
#SBATCH --mem-per-cpu=5G	# Memory per cpu in G (see also --mem-per-cpu)


# Define paths and variables
path_indexed_genome="/home/ggomez/tfm_paula/aditional_info/index/genome_index_hisat2"		# Path to HISAT2 genome index
reference="/home/ggomez/tfm_paula/librerias/archivos_TFM_Paula/ChineseLong_genome_v3.fa"	# Reference genome for index
input_dir="/home/ggomez/tfm_paula/results/bowtie_filtrado/filtered_pairs"			# Directory with .fq files
output_dir_al="/home/ggomez/tfm_paula/results/hisat2/aligned"					# Output directory for aligned files
output_dir_un="/home/ggomez/tfm_paula/results/hisat2/unaligned"					# Output directory for unaligned files
sam_outputdir="/home/ggomez/tfm_paula/results/hisat2/sam_output"				# Output directory for SAM files
bam_outputdir="/home/ggomez/tfm_paula/results/hisat2/bam_output"				# Output directory for BAM files
num_threads=16											# Number of threads for HISAT2

# Create directories if they don't exist
mkdir -p "$output_dir_al"
mkdir -p "$output_dir_un"
mkdir -p "$path_indexed_genome"
mkdir -p "$sam_outputdir"
mkdir -p "$bam_outputdir"

# Check if the index file already exists
if [[ -f "$path_indexed_genome/genome_index_hisat2.1.ht2" ]]; then
    echo "The index already exists at: $path_indexed_genome/genome_index_hisat2"
    echo "Continuing with existing index..."
else
    echo "Index does not exist. Creating genome index..."
    hisat2-build "$reference" "$path_indexed_genome/genome_index_hisat2"
    echo "Indexed genome saved at: $path_indexed_genome"
fi

# Iterate over the R1 files in the input directory
for r1_file in "$input_dir"/*_1_paired_filt.fq; do
    if [[ -f "$r1_file" ]]; then
        # Get the base name of the R1 file
        base_name=$(basename "$r1_file" _1_paired_filt.fq)

        # Build the path for the corresponding R2 file
        r2_file="$input_dir/${base_name}_2_paired_filt.fq"
        
        # Check if both R1 and R2 files exist
        if [[ -f "$r2_file" ]]; then
            # Run the hisat2 command to perform the alignment
            srun -N 1 -n 1 -c$SLURM_CPUS_PER_TASK -Q --exclusive hisat2 -x "$path_indexed_genome/genome_index_hisat2" \
                -1 "$r1_file" \
                -2 "$r2_file" \
                -S "$sam_outputdir/$base_name.sam" \
                --un-conc "$output_dir_un/$base_name" \
                --al-conc "$output_dir_al/$base_name" \
                -p "$num_threads" &
                
            echo "Processed files R1: $r1_file and R2: $r2_file"
            echo "Results saved at: $sam_outputdir/$base_name.sam"
            echo
        fi
    fi
done
wait

# Rename the outputs from --un-conc and --al-conc to match the expected extension
for r1_file in "$input_dir"/*_1_paired_filt.fq; do
    if [[ -f "$r1_file" ]]; then
        # Get the base name of the R1 file
        base_name=$(basename "$r1_file" _1_paired_filt.fq)

        # Change the extension of --un-conc and --al-conc outputs to ".fq"
        mv "$output_dir_al/${base_name}.1" "$output_dir_al/${base_name}_1_paired_filt.fq"
        mv "$output_dir_al/${base_name}.2" "$output_dir_al/${base_name}_2_paired_filt.fq"
        mv "$output_dir_un/${base_name}.1" "$output_dir_un/${base_name}_1_paired_filt.fq"
        mv "$output_dir_un/${base_name}.2" "$output_dir_un/${base_name}_2_paired_filt.fq"

        echo "Renamed output files for $base_name"
    fi
done

# Convert SAM files to BAM, sort and index them
for sam_file in "$sam_outputdir"/*.sam; do
    if [[ -f "$sam_file" ]]; then
        # Get the base name of the SAM file
        base_name=$(basename "$sam_file" .sam)
        
        # Convert SAM to BAM
        bam_file="$bam_outputdir/$base_name.bam"
        sorted_bam_file="$bam_outputdir/${base_name}_sorted.bam"
        
        echo "Converting $sam_file to $bam_file..."
        samtools view -@ $num_threads -b "$sam_file" -o "$bam_file"
        
        # Sort BAM
        echo "Sorting $bam_file to $sorted_bam_file..."
        samtools sort -@ $num_threads -o "$sorted_bam_file" "$bam_file"
        
        # Index sorted BAM
        echo "Indexing $sorted_bam_file..."
        samtools index -@ $num_threads "$sorted_bam_file"
        
        echo "Completed processing for $sam_file. Sorted and indexed BAM available at $bam_outputdir"
    fi
done

echo "Process completed. Aligned and processed files have been saved at $bam_outputdir."
exit 0

