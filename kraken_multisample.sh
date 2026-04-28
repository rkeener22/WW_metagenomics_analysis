for sample in /dodo/rk167/metafastp/norRNA/*.fastq.gz; do
	    # Get just the filename (e.g., sample1.fastq.gz) instead of the full path
	base=$(basename "$sample")
    
	kraken2 --db /home/dbs/Kraken2/k2_core_nt \
                --memory-mapping \
		--threads 30 \
		--use-names \
		--report-minimizer-data \
		--output "/dodo/rk167/krakenreads/DDnorRNA_CL0/kraken_${base}.txt" \
		--report "/dodo/rk167/krakenreports_CL0/krakenreport_${base}.txt" \
		"$sample"
done
