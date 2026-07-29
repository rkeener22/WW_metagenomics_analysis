for sample in /dodo/rk167/metafastp/norRNA/*_dedup_noribo.fastq.gz;

do

s=$(echo $sample | sed -E "s/\_.*//");

megahit -r /dodo/rk167/metafastp/norRNA/"$s"_dedup_noribo.fastq.gz --presets meta-large -o /dodo/rk167/rk_contigs/"$sample"_contigs.fasta/; 
done
