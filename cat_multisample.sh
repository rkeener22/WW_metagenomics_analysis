#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Create an output directory to keep raw data clean
mkdir -p concatenated_out

# Loop through unique sample prefixes by targeting Read 1 files
for r1_file in *_L001_R1.fastq.gz; do
    
    # Extract the unique Sample ID base name
    # Example: transforms "SampleA_L001_R1.fastq.gz" into "SampleA"
    sample_id=$(basename "$r1_file" _L001_R1.fastq.gz)
    
    echo "Processing sample: ${sample_id}"
    
    # 1. Concatenate Read 1 files in strict alphabetical/lane order
    # The wildcard '*' naturally expands in sorted alphanumeric order (L001, L002...)
    cat "${sample_id}"_L00*_R1.fastq.gz > "concatenated_out/${sample_id}_R1.fastq.gz"
    
    # 2. Concatenate Read 2 files in the exact same order
    cat "${sample_id}"_L00*_R2.fastq.gz > "concatenated_out/${sample_id}_R2.fastq.gz"

done

echo "All samples successfully concatenated!"
