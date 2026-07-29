for sample in /dodo/rk167/metafastp/norRNA/*_dedup_noribo.fastq.gz;

do

f=$(echo $sample | sed -E "s/\_.*//");

megahit -1 "$f"_1.fastq -2 "$f"_2.fastq -o "$file"_assembled.fasta/; 
done
