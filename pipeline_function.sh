#i want to make it so I can give sample name inputs and it will run!
#perhaps giving it a file header like TX_Houston_69th_20260615_S21_L001_R1_001.fastq.gz minus everything after the lane and then it will do everything from start to finish
#add dependencies in notes
#add some options?


rk_metaG() {
    local SAMPLE_ID=$1
    local R1=$2
    local R2=$3
    local OUTDIR="results_${SAMPLE_ID}"

    # Create output directory
    mkdir -p "$OUTDIR"

    # Step 1: concatenate, make it optional

    # Step 2: Quality trimming with Trimmomatic ----- make fastp
    echo "Processing QC for $SAMPLE_ID..."
    trimmomatic PE -threads 4 "$R1" "$R2" \
        "$OUTDIR/${SAMPLE_ID}_R1_paired.fastq.gz" "$OUTDIR/${SAMPLE_ID}_R1_unpaired.fastq.gz" \
        "$OUTDIR/${SAMPLE_ID}_R2_paired.fastq.gz" "$OUTDIR/${SAMPLE_ID}_R2_unpaired.fastq.gz" \
        SLIDINGWINDOW:4:20 MINLEN:50

    # Step 3: ribodetector

    # Step 4: Taxonomic classification with Kraken2
    echo "Running taxonomic classification..."
    kraken2 -db /path/to/kraken2_db \
        --paired "$OUTDIR/${SAMPLE_ID}_R1_paired.fastq.gz" \
        "$OUTDIR/${SAMPLE_ID}_R2_paired.fastq.gz" \
        --report "$OUTDIR/${SAMPLE_ID}.report.txt" \
        --output "$OUTDIR/${SAMPLE_ID}.kraken.out"

    # Maybe add in bracken here
    # Step 5: BWA test against genome file

    echo "Pipeline finished for $SAMPLE_ID!"
}
