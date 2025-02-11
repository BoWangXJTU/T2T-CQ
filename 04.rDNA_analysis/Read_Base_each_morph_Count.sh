#!/bin/bash

# Set BAM file name
bam="mat/pat_ont_sort.bam"

# Output header for the statistics table
echo -e "Chromosome\tRead_Count\tBase_Count\tFile_Name" > result_table.tsv

# Define chromosome and region information
declare -A regions=(
    ["chr13"]="1-46970"
    ["chr14"]="1-44272"
    ["chr15a"]="1-42527"
    ["chr15b"]="1-44385"
    ["chr15c"]="1-44929"
    ["chr21a"]="1-44787"
    ["chr21b"]="1-44385"
    ["chr22"]="1-43376"
)

# Iterate through each chromosome, extract sequences, and generate FASTQ files
for chr in "${!regions[@]}"; do
    region="${regions[$chr]}"
    fastq_file="$(basename $bam _sort.bam)_${chr}_sort.fastq"
    
    echo "Processing $chr ($region)..."
    samtools view "$bam" "$chr:$region" -b | samtools fastq - > "$fastq_file"
    
    if [[ -f $fastq_file ]]; then
        # Use seqtk to calculate Read and Base counts
        stats=$(seqtk size "$fastq_file" 2>/dev/null)
        if [[ -n "$stats" ]]; then
            read_count=$(echo "$stats" | awk '{print $1}')  # First column is Read count
            base_count=$(echo "$stats" | awk '{print $2}')  # Second column is Base count
            echo -e "${chr}\t${read_count}\t${base_count}\t${fastq_file}" >> result_table.tsv
        else
            echo -e "${chr}\t0\t0\t${fastq_file} (Error in seqtk size)" >> result_table.tsv
        fi
    else
        echo -e "${chr}\t0\t0\t${fastq_file} (Not Found)" >> result_table.tsv
    fi
done

echo "All results have been saved to result_table.tsv"
