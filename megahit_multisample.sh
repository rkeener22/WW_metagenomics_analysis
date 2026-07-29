for sample in /dodo/rk167/metafastp/norRNA/*_dedup_noribo.fastq.gz;

do

s=$(echo $sample | sed -E "s/\_.*//");

megahit -r "$s"_dedup_noribo.fastq.gz -o "$sample"_assembled.fasta/; 
done
