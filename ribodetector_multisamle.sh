for file in input_dir/*.fastq.gz; do
    base=$(basename "$file" .fastq.gz)
    ribodetector_cpu -t 10 -l 100 \
        -i "$file" \
        -o "output_dir/${base}_non_rrna.fastq.gz" \
        -e rrna
done
