###  RNA-Seq: Preparation and FASTQ File Download Guide

> **Before You Begin**  
> Ensure you are using the **same data source** consistently throughout your RNA-Seq analysis to avoid inconsistencies in results.

---------

###  Downloading FASTQ Files from ENA Using `wget`

If you want to download data from the **European Nucleotide Archive (ENA)** using `wget`, follow these steps:

###  Step-by-Step Guide

1. Visit the ENA browser: [https://www.ebi.ac.uk/ena/browser/home](https://www.ebi.ac.uk/ena/browser/home)
2. Search for your **SRR accession** (e.g., `SRR1551114`)
3. Scroll down to the **"FASTQ Files"** section
4. Right-click and copy the download links for:
   - `*_1.fastq.gz` (Read 1)
   - `*_2.fastq.gz` (Read 2)
5. Use `wget` to download the files:

```bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR155/004/SRR1551114/SRR1551114_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR155/004/SRR1551114/SRR1551114_2.fastq.gz
```

------------------

### Unzip and save the uncompressed file

This will create a plain .fastq file (warning: could be huge).
gunzip SRR33228995.fastq.gz

If you want to keep both compressed and uncompressed files: gzip -c (source directory) > (destination directory)
```
gunzip -c ~/datasets/SRR1551114/zipped/SRR1551114_2.fastq.gz > ~/datasets/SRR1551114/raw/SRR1551114_2.fastq
```
- -c : you must provide a single input file.
- '>' : keep original .gz file

------------------------------------------

### RNA-Seq: Downloading FASTQ Files from NCBI


If you want to download from NCBI, Use SRA ( Sequence Read Archieve ) -Tool

Tool: sra-toolkit
check command: vdb-dump --version

e.g., 
```
fasterq-dump SRR33542395 --split-files -O ~/Downloads --temp ~/Downloads/tmp --progress
```
- --split-siles : downloads both reads if available
-  -o: destination for saving the file
-  --temp ~/Downloads/tml: temporary directory to record the progress as too mucch large data could look like stuck while its downloading in backend
-  --progress : to visualize the downloading progress





------------------------------------------------------------------------------------
#   🧬 STARTING RNA-SEQ ANALYSIS 🧬                   
------------------------------------------------------------------------------------

## ▶️ Run FastQC on the downloaded files:

```
fastqc SRR1551114_1.fastq.gz SRR1551114_2.fastq.gz -o ~/fastq_qc_reports
```
 or
 ```
fastqc -t 4 -o /home/vzscyborg/rnaseq/fastqc_rep /home/vzcyborg/datasets/rnads/SRR1551114_1.fastq.gz
```
- fastqc - command calling the quality check function
- -t 4 - t refers to the threads using depanding upon the per core threads
- -o refers the path where the report of the quality check will save followed by the proper path
- SRRxxxxxxxxx.fastq.gz - files on which quality check will run

--------------------------------------

## ▶️ Run Fastp, if the Fastqc report reaches quality failure

Fastp require to ]run if following conditions arises:
  1. Adapter contamination 
  2.  Low-quality tiles 
  3. Biased base content 
  4. verrepresented sequences 

```
fastp -i /home/vzscyborg/ngs/datasets/SRR33542395 -o /home/vzscyborg/ngs/output/SRR33542395/fp/SRR33542395_clean.fastq --detect_adapter_for_pe --html /home/vzscyborg/ngs/output/SRR33542395/fp/SRR33542395_clean.html --json /home/vzscyborg/ngs/output/SRR33542395/fp/SRR33542395_clean.json --thread 1
```
- -i: for input file followed by file directory, if two read second read start with \-I (source directory) 
-  -o: for output file followed by file directory with extension, if two out reads second read start with /-O (destination directory) 
-  --html : for making the output report in html formate along with the file name 
-  --thread define the number of thread fastp should use

----------------------------------------------------------

## ▶️ Run hisat2

hisat2 have two functions, one is for making the genomic index of the reference genome.
Another function is to do refernce genomce alignment

### 1. Genome Indexing: 
 For reference genome i took GRCm39 from NCBI ( Genome ), downloaded the Genome Sequences (FASTA)
 Following comnand is for the generating the genemic index
```
  hisat2-build /home/vzscyborg/ngs/datasets/gtf/GRCm39.primary_assembly.genome.fa /home/vzscyborg/ngs/output/SRR33542395/hisat/indx/GRCm39_index
```
-  hisat2-build - 	The HISAT2 utility that builds the index files from a reference genome (FASTA format) 
-  /home/vzscyborb... - path to the reference genome 
-  grcm_index - This is the base name for the output index files. HISAT2 will generate 8 files named grcm_index.1.ht2, grcm_index.2.ht2, ..., grcm_index.8.ht2.

### 2. ALignment 
 For the alignment require multiple files for 2 input method one for genome index and another is cleaned fastq file generated from fastp
```
hisat2 -p 1 -x /home/vzscyborg/ngs/output/SRR33542395/hisat/indx/grcm_index -U /home/vzscyborg/ngs/output/SRR33542395/fp/SRR33542395_clean.fastq -S /home/vzscyborg/ngs/output/SRR33542395/hisat/align/SRR33542395_hisat.sam
```
-  hisat2 : HISAT2 is a fast and sensitive alignment tool for mapping sequencing reads (usually RNA-seq) to a reference genome. 
- -x : Specifies the basename of the HISAT2 index for the reference genome to which you want to align your reads. This is the index you built earlier from the mouse genome (GRCm39). 
-  -p 1 : Use 1 CPU thread for the alignment process. This limits HISAT2 to run on a single processor core 
-  -U : Specifies the input file containing unpaired single-end reads in FASTQ format. This is the cleaned sequencing reads you want to align. 
-  -S : Specifies the output file where the aligner will write the results in SAM format (Sequence Alignment/Map). This file contains detailed information on how each read aligns to the reference genome

------------------------------------


## ▶️ Run samtools

### 1. convers the SRR35542395.sam to SRR35542395.bam
```
   samtools view -@ 3 -s 0.5 -b /home/vzscyborg/ngs/output/SRR33542395/hisat/align/SRR33542395_hisat.sam > /home/vzscyborg/ngs/output/SRR33542395/smtl/SRR33542395.bam
```
- samtools view : Converts between different formats (SAM ⇄ BAM), filters, and subsamples reads 
- -@ : 	Use 3 CPU threads to speed up the conversion. 
- -s 0.5 : 	Subsample 50% of reads randomly. The 0.5 means each read has a 50% chance of being retained. This is useful for downsampling large datasets  ( OPTIONAL )
- -b Output in BAM format (binary, compressed form of SAM). 
- 	Redirects output to a file.

### 2. BAM file soritng
```
  samtools sort -@ 3 -o /home/vzscyborg/ngs/output/SRR33542395/smtl/SRR35542395_sorted.bam /home/vzscyborg/ngs/output/SRR33542395/smtl/SRR33542395.bam
```
-  samtools sort : Convert Bam file into sorted bam file.
-  -@ : Use 3 CPU threads to speed up the conversion.
-  -o : specify the output file followed by the output file path along with the filename with .bam file extension
-  in last the directiory to the input bam file

### 3. Indexing of sorted BAM file
```
samtools index -@ 3 -o /home/vzscyborg/ngs/output/SRR33542395/smtl/SRR33542395.bai /home/vzscyborg/ngs/output/SRR33542395/smtl/SRR33542395_sorted.bam
```
-  samtools index : Convert Sorted  Bam file into indexed bam file.
-  -@ : Use 3 CPU threads to speed up the conversion.
-  -o : specify the output file followed by the output file path along with the filename with .bam file extension
-  in last the directiory to the input sorted bam file


### or all three steps in one time:-
```
samtools view -@ 1 -bS /home/vzscyborg/ngs/output/SRR33542395/hisat/align/SRR33542395_hisat.sam | \
samtools sort -@ 1 -o /home/vzscyborg/ngs/output/SRR33542395/smtl/SRR33542395_sorted.bam
samtools index -@ 1 /home/vzscyborg/ngs/output/SRR33542395/smtl/SRR33542395_sorted.bam
```

------------------------------------

## ▶️ Run FeatureCounts

for creating the count matrix of the gene/exons must have genemoe assembaly downloaded from GENCODE ( recommanded ) and dataset set must matches the reference genome.
