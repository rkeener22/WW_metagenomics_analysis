# Test metagenomics SOP

## 1. Getting things ready!

### Accessing the Stadler lab section of the p00 Server
first you will need a linux/unix enabled "portal" to your computers terminal. For apple users this is already available but for windows users, consider Ubuntu or Gitbash. There are also more built out options that are more user friendly, such as the terminal tab in r studio. 

Open the terminal and type ```ssh YOURRICEID@p00.cs.rice.edu```. You will be prompted for your normal rice password. Once in, you will be located in your home directory, usually ```home/users/YOURUSERID```. We also have a data/ and dodo/ directory. we will do most of our work on dodo/ due to space constraints on the other drives. 

### Creating a Conda Environment & Adding needed programs
BWA, Kraken2, krakentools, fastp, ribodetector, ncbi_datasets, megahit

```conda create -n metaG -y bwa kraken2 krakentools fastp ribodetector ncbi-datasets-cli megahit```
or mamba
```mamba create -n metaG -y bwa kraken2 krakentools fastp ribodetector ncbi-datasets-cli megahit```

activate environment anywhere, including in a tmux session so you can close the terminal!

```conda activate metaG```


## 2. Download data - we won't do this today

All of the raw data is saved in ```/dodo/stadler_ww/RAWFILES``` under the date of the run and project code given by O'Connor lab. To download:

### From Missouri sequencing core

While in the p00 server in the ```/dodo/stadler_ww/RAWFILES``` folder (or folder of your choice), enter the following command and then the password when prompted. This is best done in a tmux terminal multiplier so you can leave it and it wont impact download. 

```sftp ext-ls58@download.bioinformatics.missouri.edu```

get command will download whatever needed to the location you started in p00, add -r for whole folders

```get -r folder_of_choice/```




## 3. Concatenate 

Raw files are located in the ```/dodo/stadler_ww/SOP_test/``` folder. Navigate here with ```cd /dodo/stadler_ww/SOP_test/```

Next you will use ```cat``` to create a single file for each paired end. This concatenates (combines end-to-end) all files (usually 8 lanes for each date, ive shrunk it to two for this test) from each paired end (R1 & R2). You must ensure these are in order so that the R1 and R2 files will be consistent for later processing. Generally we will want to keep the header of the file the same when creating new files, everything that comes before the lane number. 

R1
```Cat 001_R1.fastq.gz 002_R1.fastq.gz > combined_R1.fastq.gz```
R2
```Cat 001_R1.fastq.gz 002_R1.fastq.gz > combined_R1.fastq.gz```



## 4. Data cleaning
### Quality filtering
Here, we use ```fastp``` to remove low quality reads, reads with too many Ns, remove reads shorter than 100 bp, merge paired ends for longer sequence and reduced file size, repair low quality and N bases with other paired end, deduplicate PCR replicates. Here, there are several possibilities for cleaning the data. 

#### My run
For initial merging and cleaning as well as statistics:
```fastp -m -i combined_R1.fastq.gz -I combined_R2.fastq.gz --merged_out combined.fastq.gz --dedup -l 100 -w 10```




## 5. Ribosomal RNA removal

Due to the rna-biased extraction method, bacterial rRNA is highly present in the data. These reads are difficult to assign to specific species, so they are removed for downstream assignment. These reads may be important in genome assembly or other tools. Here, i used a database-free neural network deep learning model called ```ribodetector``` to remove known and possible novel rRNA sequences from each sequencing file. This is faster on GPU but can be run on CPU with many threads. **100,000,000 reads takes about 1 hour on 40 threads**.  

```ribodetector_cpu -l 150 -i combined.fastq.gz -o _dedup_noribo_combined.fastq.gz -t 20 --chunk_size 2000```



## 6. Taxonomic Classification (and Human Read Assignment)
### Using kmer-based classifier for efficient assignment
Here, I use ```Kraken2``` to quickly assign reads to a core_nt dataset of refseq sequences. It includes all taxonomies including fungi, viruses, bacteria, and eukaryotes. I find that using a confidence level of 0.1 has been effective at stringency but also still getting to the species level much of the time. This produces a report which lists the percent, raw read count, and unique identifiers assigned to each taxa as well as a list of the reads and their assignment. 

```kraken2 --db /home/dbs/Kraken2/k2_core_nt --threads 10 _dedup_noribo_combined.fastq.gz --confidence 0.1 --report krakenreport_01.txt --report-minimizer-data --use-names > krakenreads_01.txt```

### Combine reports for data analysis - we don't need to do this today

Using ```krakentools```, we will manipulate the files for downstream analysis.

This provides a single file with level and total raw read numbers for each date and taxa. 

```combine_kreports.py -r *dedup_noribo_C01.txt -o Combined.txt```

Following this, you can download and clean this file for downstream analysis. This will include removing the header, dividing by total reads per sample and multiplying that by 1 million, and possibly normalizing to flow rate data. 

#### Extracting reads for alignments - we don't need to do this today

This extracts reads from specific taxa for alignment or reports. 

```extract_kraken_reads.py -k krakenreads/krakenreads_combined.txt -s metafastp/combined.fastq.gz -o bpertreads.fasta -t 520```

```-k``` : kraken-assigned read file for sample of interest. 

```-s``` : Original merged read fastq file.

```-o``` : Read output file. 

```-t``` : list of taxa to extract. 


#### Download or upload from your computer
we will download our file to your computer to look through it in excel! I load it to my downloads folder but you can do another place if you feel like it. 

```scp USERID@p00.cs.rice.edu:/dodo/stadler_ww/SOP_test/yourkrakenreport_file.txt ~/Downloads/```



# We won't do this today

## 7. Align sequences to genome of interest - we wont do this today
Using the reads obtained above, we will utilize alignment platforms and blast to confirm the reads. 

### Using BWA 

#### preparing genomes
First, you will need to download refseq genomes of your species of interest from NCBI.

```conda activate ncbi_datasets```
```datasets download genome taxon "Bordetella" --reference --filename /dodo/rk167/bordatella_ref.zip```
```unzip bordatella_ref.zip -d Bordetella_ref```

Next, you will need to concatenate the genomes if they arent already using ```cat```.

Finally, you will need to index the genomes. 

```./bwa/bwa index ref.fna```

#### genome alignment
Next, you will take the .fasta file created from the read extraction above and align it to the genomes of interest.

```./bwa/bwa mem /dodo/rk167/Bordetella_ref/Bordetella_combined.fna reads.fastq.gz > reads_aligned.sam```

Then you will need to turn your .sam file into a sorted .bam

```samtools view -bS input.sam | samtools sort -o sorted_output.bam –```

You can check how they look 

```samtools view -F 0x4 -q 1 your_file.sorted.bam```

then index it

```samtools index ```

and upload your indexed genome and reads into IGV. 


samtools view -S -b output.sam > output.bam
samtools sort output.bam -o sorted.bam
samtools index sorted.bam


