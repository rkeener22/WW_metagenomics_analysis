# 1. Download data

All of the raw data is saved in ```/dodo/stadler_ww/date_of_sequencing```

## From Missouri sequencing core

While in the p00 server in the ```/dodo/stadler_ww``` folder (or folder of your choice), enter the following command and then the password given to you when prompted. This is best done in a tmux terminal multiplier so you can leave it and it wont impact download. 

```sftp ext-ls58@download.bioinformatics.missouri.edu```

get command will download whatever needed to the location you started in p00, add -r for whole folders

```get -r folder_of_choice/```

## From Globus (not recommended but sometimes necessary)

This may be necessary if downloading data from another lab (no longer on sequencing core server) using the **Globus webapp** logged into your rice account or whatever account the sharepoint was added to. From secure endpoint **SecureBio-Houston Data** to compute cluster **Rice University CRC RDF Collection** Globus endpoint (under ```ls57``` folder, you will need permission to use this account and it can only hold 500gb at a time) then to globusconnectpersonal endpoint on p00 server (you will have to set this up yourself, easy directions online). Use the following code to turn on personal endpoint on p00.

```/globusconnectpersonal -start &```

You can then add files to the server using the **Globus webapp**

# 2. Concatenate 
Make a new directory in ```/data/stadler_ww``` with the same header as the raw data folder.

```mkdir /dodo/stadler_ww/date_of_sequencing```

Next you will use ```cat``` to create a single file for each paired end. This concatenates (combines end-to-end) all files (usually 8 lanes for each date) from each paired end (R1 & R2). You must ensure these are in order so that the R1 and R2 files will be consistent for later processing.

```Cat 001_R1.fastq.gz 002_R1.fastq.gz 003_R1.fastq.gz 004_R1.fastq.gz 005_R1.fastq.gz 006_R1.fastq.gz 007_R1.fastq.gz 008_R1.fastq.gz > /data/stadler_ww/date_of_sequencing/combined_R1.fastq.gz```

 Find cat files in ```/data/stadler_ww/date_of_sequencing```.

# 3. Data cleaning
## Quality filtering
Here, we use ```fastp``` to remove low quality reads, reads with too many Ns, remove reads shorter than 100 bp, merge paired ends for longer sequence and reduced file size, repair low quality and N bases with other paired end, deduplicate PCR replicates. Here, there are several possibilities for cleaning the data. 

### My run
For initial merging and cleaning as well as statistics:
```fastp -m -i combined_R1.fastq.gz -I combined_R2.fastq.gz --merged_out combined.fastq.gz -l 100 -w 10```

cleaned and merged fastp files are found in ```/dodo/rk167/metafastp```

For deduplication: 
```fastp -i input.fastq.gz -o dedupe.fastq.gz -D -w 10```
or start with
```fastp -m -i combined_R1.fastq.gz -I combined_R2.fastq.gz --merged_out combined.fastq.gz --dedup -l 100 -w 10```

Deduplicated files are saved in ```/dodo/rk167/metafastp/dedupe```


#### Defaults include: 

Adapter Trimming: Detects common adapters automatically, including Illumina and Nextera.

Quality Trimming: Trims low-quality bases from both ends (sliding window from 5' to 3') with a default threshold of Q15. 40% (int [=40]) of bases unqualified is max. if one read's number of N base is >5, then this read/pair is discarded.

Length Filtering: Reads shorter than 15bp are discarded.

Base Correction: Enabled for paired-end data. If a good overlap exists, mismatching bases are corrected if one has a high-quality score (>Q30) and the other is low-quality (<Q15).

Threading: Uses 4 threads by default.

Duplication Evaluation: fastp evaluates duplication rate, and this module may use 1G memory and take 10% ~ 20% more running time.

#### Variables I set include:

```-m``` : merge reads

```-i``` : input R1

```-I``` : input R2

```-l``` : The minimum length requirement (100bp here), should remove short reads **after merging**.

```-D, --dedup``` : Deduplicates the data by hash algorithm and only counts exact matches. The ```dup_calc_accuracy level``` is default to 3 wich has a hash buffer number of 3. 

```--merged_out``` : given to specify the file to store merged reads. The merged reads are also filtered.

```-w``` : changed to 8 or 10 threads but can be altered as needed

##### Note on merging reads

This may not be best practice for those methods that use genome alignment and contigs for species confirmation. 

##### Note on deduplication

```fastp``` uses a hash algorithm to find the identical sequences. Due to the possible hash collision, about 0.01% of the total reads may be wrongly recognized as deduplicated reads. Normally this may not impact the downstream analysis. The accuracy of calculating duplication can be improved by increasing the hash buffer number or enlarge the buffer size. Please refer to the deduplication table: https://github.com/opengene/fastp?tab=readme-ov-file#simple-usage


#### To include in future?
```--out1``` and ```--out2``` : will be the reads that cannot be merged successfully, but both pass all the filters.

```--include_unmerged``` : can be enabled to make reads of --out1, --out2, --unpaired1 and --unpaired2 redirected to --merged_out. So you will get a single output file. This option is disabled by default.

```-j```, ```--json``` : the json format report file name (string [=fastp.json])

## Ribosomal RNA removal

Due to the rna-biased extraction method, bacterial rRNA is highly present in the data. These reads are difficult to assign to specific species, so they are removed for downstream assignment. These reads may be important in genome assembly or other tools. Here, i used a database-free neural network deep learning model called ```ribodetector``` to remove known and possible novel rRNA sequences from each sequencing file. This is faster on GPU but can be run on CPU with many threads. **100,000 reads takes about 1 hour on 40 threads**.  

```ribodetector_cpu -l 150 -i dedupe/merged_dedup.fastq.gz -o norRNA/merged_dedup_noribo.fastq.gz -t 40 --chunk_size 2000```

The rRNA depleted data is stored in ```/dodo/rk167/metafastp/noribo```

#### Default parameters: 
Model (--model-file): model_len70_101 (packaged model)

Input/Output: Supports FASTQ/FASTA, paired or single-end.

Discordant Pair Behavior: If using none for prediction labels, discordant read pairs are discarded.

Model Accuracy: Accuracy reduces for reads shorter than 40bp.

Memory Usage: The --chunk_size parameter should be adjusted if you have low memory.

#### My added parameters:

```-l``` : Mean read length of the sample. Required, helps estimate time. 

```-i``` : Input

```-o``` : Output of only non-rRNA sequences.

```-t``` : Number of threads, I use 20 or 40. Memory will be about 50 gb for 20 threads. 

```--chunk_size``` : Number of reads to divide amoungst the threads each round. Should be a multiple of the threads. Too few or too many will increase processing time. I use 1000/20 threads. 

## Human Read Depletion
not done yet

# 4. Taxonomic Classification
## Using kmer-based classifier for efficient assignment
Here, I use ```Kraken2``` to quickly assign reads to a core_nt dataset of refseq sequences. It includes all taxonomies including fungi, viruses, bacteria, and eukaryotes. I find that using a confidence level of 0.1 has been effective at stringency but also still getting to the species level much of the time. This produces a report which lists the percent, raw read count, and unique identifiers assigned to each taxa as well as a list of the reads and their assignment. 

```kraken2 --db /home/dbs/Kraken2/k2_core_nt --threads 10 /dodo/rk167/metafastp/combined.fastq.gz --confidence 0.1 --report /dodo/rk167/krakenreports_C01/krakenreport_combined.txt --report-minimizer-data --use-names > /dodo/rk167/krakenreads/kraken_combined.txt```

reports are saved to ```/dodo/rk167/krakenreports_C01```
reads are saved to ```/dodo/rk167/krakenreads```

#### Default parameters: 
Confidence (--confidence): 0 (no confidence threshold; a single k-mer match can classify a read). I changed this.

Minimum Hit Groups (--minimum-hit-groups): 2 (minimum number of hits required to classify a read).

K-mer length (-k): 31.

Minimizer length (-l): 15.

Spaced Seed (-s): 7 (used for nucleotide databases).

Memory: Loads database into RAM (no memory mapping).

#### My additional parameters:
```--db``` : database directory. The one I used was assmbled by someone else on **XXXX**

```--threads``` : How many threads to use, default is 1 but 10 is way faster. 

```--confidence``` : Confidence level, proportion of kmers that have to map to that taxa to be assigned. I did 0.1 but also tried 0.51 (too stringent). 

```--report``` : if you want a report this is the output. 

```--report-minimizer-data``` : this provides the number of unique kmers asigned to the taxa. Useful for figuring out the coverage. High read low minimizer count? Probably not real. 

```--use-names``` : Provides scientific name of taxa in addition to the taxa number.

## Data cleaning for Analysis

Using ```krakentools```, we will manipulate the files for downstream analysis.

### Combine reports for data analysis

This provides a single file with level and total raw read numbers for each date and taxa. 

```combine_kreports.py -r *dedup_noribo_C01.txt -o Combined.txt```

Following this, you can download and clean this file for downstream analysis. This will include removing the header, dividing by total reads per sample and multiplying that by 1 million, and possibly normalizing to flow rate data. 

### Extracting reads for alignments

This extracts reads from specific taxa for alignment or reports. 

```extract_kraken_reads.py -k krakenreads/krakenreads_combined.txt -s metafastp/combined.fastq.gz -o bpertreads.fasta -t 520```

```-k``` : kraken-assigned read file for sample of interest. 

```-s``` : Original merged read fastq file.

```-o``` : Read output file. 

```-t``` : list of taxa to extract. 

#### Other settings:

```--exclude``` : Instead of finding reads matching specified taxids, finds reads NOT matching specified taxids.

```-r, --report MYFILE.KREPORT``` : Kraken report file (required if specifying --include-children or --include-parents)

```--include-children``` : include reads classified at more specific levels than specified taxonomy ID levels.

```--include-parents``` : include reads classified at all taxonomy levels between root and the specified taxonomy ID levels.

```--max #``` : maximum number of reads to save.

```--append``` : if output file exists, appends reads

# 5. Align sequences to genome of interest
Using the reads obtained above, we will utilize alignment platforms and blast to confirm the reads. 

## Using BWA 

First, you will need to download refseq genomes of your species of interest from NCBI. Next, you will need to concatenate the genomes if they arent already. Finally, you will need to index the genomes. 

```./bwa/bwa genomes XXXX```

Next, you will take the .fasta file created from the read extraction above and align it to the genomes of interest.

```./bwa/bwa mem hpvgenomes.fa 69street-08052025.fastq.gz > HPV082525.sam```

Then you will need to turn your .sam file into a sorted .bam

```samtools view -bS input.sam | samtools sort -o sorted_output.bam –```

You can check how they look 

```samtools view -F 0x4 -q 1 your_file.sorted.bam```

then index it

```samtools index ```

and upload your indexed genome and reads into IGV. 

