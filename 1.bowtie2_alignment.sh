#!/bin/bash

#SBATCH --job-name=GetBowtieReads      # Job name to display with squeue
#SBATCH --output=FileJob_%j.out        # Output file
#SBATCH --ntasks=18                    # Number of tasks
#SBATCH --cpus-per-task=2              # Number of CPUs per task
#SBATCH -t 1-00:00:00                  # Runtime in days-hours:minutes:seconds
#SBATCH --qos short                    # The QoS to submit the job
#SBATCH --mem-per-cpu=5G               # Memory per CPU in GB (see also --mem-per-cpu)

input_dir="/home/ggomez/tfm_paula/librerias/archivos_TFM_Paula"
output_dir_al="/home/ggomez/tfm_paula/results/bowtie/aligned"
output_dir_un="/home/ggomez/tfm_paula/results/bowtie/unaligned"
path_indexed_genome="/home/ggomez/tfm_paula/aditional_info/index"
reference="/home/ggomez/tfm_paula/aditional_info/back_splicing_sequences.fasta"
sam_outputdir="/home/ggomez/tfm_paula/results/bowtie/sam_output"
num_threads=18

mkdir -p "$output_dir_al"
mkdir -p "$output_dir_un"
mkdir -p "$path_indexed_genome"
mkdir -p "$sam_outputdir"

module load bowtie2

# Index the back splicing sequences
bowtie2-build "$reference" "$path_indexed_genome/back_splicing_sequences"
echo "Indexed genome saved at: $path_indexed_genome"

# Iterate over R1 files in the input directory
for r1_file in "$input_dir"/*_1_paired.fq; do
    if [[ -f "$r1_file" ]]; then
        # Get the base name of the R1 file
        base_name=$(basename "$r1_file" _1_paired.fq)

        # Build the corresponding R2 file path
        r2_file="$input_dir/${base_name}_2_paired.fq"
        
        # Check if both R1 and R2 files exist
        if [[ -f "$r2_file" ]]; then
            # Run bowtie2 to perform alignment
            srun -N 1 -n 1 -c$SLURM_CPUS_PER_TASK -Q --exclusive bowtie2 -x "$path_indexed_genome/back_splicing_sequences" \
                -q -1 "$r1_file" \
                -2 "$r2_file" \
                -S "$sam_outputdir/$base_name.sam" \
                --un-conc "$output_dir_un/$base_name" \
                --al-conc "$output_dir_al/$base_name" \
                --threads "$num_threads" &
                
            echo "Processed files R1: $r1_file and R2: $r2_file"
            echo "Results saved in: $output_dir_al/$base_name.sam, and $output_dir_al/$base_name"
            echo
        fi
    fi
done
wait
          
# Iterate again to rename the output files from --un-conc and --al-conc to ".fq"
for r1_file in "$input_dir"/*_1_paired.fq; do
    if [[ -f "$r1_file" ]]; then
        # Get the base name of the R1 file
        base_name=$(basename "$r1_file" _1_paired.fq)
       
        # Rename output files to have .fq extension
        mv "$output_dir_al/${base_name}.1" "$output_dir_al/${base_name}_1_paired.fq"
        mv "$output_dir_al/${base_name}.2" "$output_dir_al/${base_name}_2_paired.fq"
        mv "$output_dir_un/${base_name}.1" "$output_dir_un/${base_name}_1_paired.fq"
        mv "$output_dir_un/${base_name}.2" "$output_dir_un/${base_name}_2_paired.fq"
       
        echo "Processed file $base_name"
    fi
done

echo "Process completed. Aligned files have been saved in $output_dir."
exit 0

