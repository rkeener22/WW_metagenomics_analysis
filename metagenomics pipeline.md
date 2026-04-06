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
## fastp
Here, we use fastp to remove low quality reads, reads with too many Ns, remove reads shorter than 100 bp, merge paired ends for longer sequence and reduced file size, repair low quality and N bases with other paired end, deduplicate PCR replicates. Here, there are several possibilities for cleaning the data. 

### My run

```fastp -m -i combined_R1.fastq.gz -I combined_R2.fastq.gz --merged_out combined.fastq.gz --dedup -l 100 -w 10```

#### defaults include: 

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
```--dedup``` : Deduplicates the data by hash algorithm and only counts exact matches. The ```dup_calc_accuracy level``` is default to 3 wich has a hash buffer number of 3. 
```--merged_out``` : given to specify the file to store merged reads. The merged reads are also filtered.
```-w``` : changed to 8 or 10 threads but can be altered as needed

##### Note on merging reads

This may not be best practice for those methods that use genome alignment and contigs for species confirmation. 

##### Note on deduplication

```fastp``` uses a hash algorithm to find the identical sequences. Due to the possible hash collision, about 0.01% of the total reads may be wrongly recognized as deduplicated reads. Normally this may not impact the downstream analysis. The accuracy of calculating duplication can be improved by increasing the hash buffer number or enlarge the buffer size. Please refer to the deduplication table: https://github.com/opengene/fastp?tab=readme-ov-file#simple-usage



#### To include in future?
--out1 and --out2 will be the reads that cannot be merged successfully, but both pass all the filters.
--include_unmerged can be enabled to make reads of --out1, --out2, --unpaired1 and --unpaired2 redirected to --merged_out. So you will get a single output file. This option is disabled by default.


