#!/bin/bash

# Create output folder if it does not exist
mkdir -p trimmed_reads

# Loop through all Forward (R1) read files
for r1_file in *_R1.fastq.gz; do
    # Dynamically find the matching Reverse (R2) file name
    r2_file="${r1_file/_R1.fastq.gz/_R2.fastq.gz}"
    
    # Extract the sample prefix name for clean reporting
    sample_name=$(basename "$r1_file" _R1.fastq.gz)
    
    echo "Processing paired-end sample: ${sample_name}"
    
    # Execute fastp
    fastp \
      -i "$r1_file" \
      -I "$r2_file" \
      -o "trimmed_reads/${sample_name}_trimmed_R1.fastq.gz" \
      -O "trimmed_reads/${sample_name}_trimmed_R2.fastq.gz" \
      -h "trimmed_reads/${sample_name}_report.html" \
      -j "trimmed_reads/${sample_name}_report.json" \
      --thread 4
done
