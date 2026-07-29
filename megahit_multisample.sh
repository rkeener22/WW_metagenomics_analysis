for file in *_1.fastq;

do

f=$(echo $file | sed -E "s/\_.*//");

megahit -1 "$f"_1.fastq -2 "$f"_2.fastq -o "$file"_assembled.fasta/; 
done
