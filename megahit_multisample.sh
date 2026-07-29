for sample in /dodo/rk167/metafastp/norRNA/*_dedup_noribo.fastq.gz; do
    # Extract just the filename (removes directory path)
    filename=$(basename "$sample")
    
    # Extract just the '*' part (removes '_dedup_noribo.fastq.gz')
    s="${filename%_dedup_noribo.fastq.gz}"

    # Run MEGAhit using the extracted sample name
    megahit -r "$sample" --presets meta-large -o "/dodo/rk167/rk_contigs/${s}_contigs"
done
