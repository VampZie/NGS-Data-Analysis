<p align="center">
  <img src="https://capsule-render.vercel.app/api?type=waving&color=0:e96443,100:904e95&height=180&section=header&text=INSTRUCTIONS&fontSize=40&fontColor=fff&animation=twinkling" alt="Waving Banner"/>
</p>

---

> **Before You Begin**  
> Always use the **same data source** throughout your RNA-Seq analysis to ensure data consistency and reproducible results.

---

## 🧬 Table of Contents

- [Downloading FASTQ Files from ENA](#downloading-fastq-files-from-ena-using-wget)
- [Unzipping FASTQ Files](#unzip-and-save-the-uncompressed-file)
- [Downloading FASTQ Files from NCBI](#rna-seq-downloading-fastq-files-from-ncbi)
- [Quality Control with FastQC](#️-run-fastqc-on-the-downloaded-files)
- [Read Trimming with Fastp](#️-run-fastp-if-the-fastqc-report-reaches-quality-failure)
- [Genome Alignment with HISAT2](#️-run-hisat2)
- [BAM Processing with Samtools](#️-run-samtools)
- [Counting Reads with FeatureCounts](#️-run-featurecounts)

---

# 🎯 RNA-Seq: Preparation and FASTQ File Download Guide

## Downloading FASTQ Files from ENA Using `wget`

If you need RNA-Seq raw data from the **European Nucleotide Archive (ENA)**:

### **Step-by-Step Guide**
1. Visit [ENA browser](https://www.ebi.ac.uk/ena/browser/home)
2. Search for your **SRR accession** (e.g., `SRR1551114`)
3. In the **"FASTQ Files"** section, copy download links for:
   - `*_1.fastq.gz` (Read 1)
   - `*_2.fastq.gz` (Read 2)
4. Use `wget` to download files:
```bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR155/004/SRR1551114/SRR1551114_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR155/004/SRR1551114/SRR1551114_2.fastq.gz
```

---

## Unzip and Save the Uncompressed File

Unzip the downloaded file (may be large):

```bash
gunzip SRR33228995.fastq.gz
```

To keep both compressed and uncompressed files:
```bash
gunzip -c ~/datasets/SRR1551114/zipped/SRR1551114_2.fastq.gz > ~/datasets/SRR1551114/raw/SRR1551114_2.fastq
```
- `-c` : Outputs to stdout, keeps the original file
- `'>'` : Saves output as a new file

---

## RNA-Seq: Downloading FASTQ Files from NCBI

If you want to download from NCBI, Use SRA ( Sequence Read Archieve ) -Tool

Tool: sra-toolkit
check command: vdb-dump --version

- Download using `fasterq-dump`:
  ```bash
  fasterq-dump SRR33542395 --split-files -O ~/Downloads --temp ~/Downloads/tmp --progress
  ```
- --split-siles : downloads both reads if available
-  -o: destination for saving the file
-  --temp ~/Downloads/tml: temporary directory to record the progress as too mucch large data could look like stuck while its downloading in backend
-  --progress : to visualize the downloading progress


------------------------------------------------------------------------------------
#   🧬 STARTING RNA-SEQ ANALYSIS 🧬                   
------------------------------------------------------------------------------------

## ▶️ Run FastQC on the Downloaded Files

**Run quality check:**
```bash
fastqc SRR1551114_1.fastq.gz SRR1551114_2.fastq.gz -O ~/fastq_qc_reports --progress
```
or
```
fastqc -t 4 -o fastqc_1.fastq fastqc_2.fastq  -O /home/vzscyborg/rnaseq/fastqc_rep /home/vzcyborg/datasets/rnads/SRR1551114_1.fastq.gz --progress
```
- fastqc - command calling the quality check function
- -t 4 - t refers to the threads using depanding upon the per core threads
- -O refers the path where the report of the quality check will save followed by the proper path
- SRRxxxxxxxxx.fastq.gz - files on which quality check will run

---

## ▶️ Run Fastp, if FastQC Report Fails

Fastp require to run if following conditions arises:
  1. Adapter contamination 
  2.  Low-quality tiles 
  3. Biased base content 
  4. verrepresented sequences
  5. 
```bash
fastp -i /home/vzscyborg/ngs/datasets/SRR33542395 \
      -o /home/vzscyborg/ngs/output/SRR33542395/fp/SRR33542395_clean.fastq \
      --detect_adapter_for_pe \
      --html /home/vzscyborg/ngs/output/SRR33542395/fp/SRR33542395_clean.html \
      --json /home/vzscyborg/ngs/output/SRR33542395/fp/SRR33542395_clean.json \
      --thread 1
```
- -i: for input file followed by file directory, if two read second read start with \-I (source directory) 
-  -o: for output file followed by file directory with extension, if two out reads second read start with /-O (destination directory) 
-  --html : for making the output report in html formate along with the file name 
-  --thread define the number of thread fastp should use


---

## ▶️ Run HISAT2

hisat2 have two functions, one is for making the genomic index of the reference genome.
Another function is to do refernce genomce alignment

### 1. **Genome Indexing**
 For reference genome i took GRCm39 from NCBI ( Genome ), downloaded the Genome Sequences (FASTA)
 Following comnand is for the generating the genemic index
 
```bash
hisat2-build /home/vzscyborg/ngs/datasets/gtf/GRCm39.primary_assembly.genome.fa \
             /home/vzscyborg/ngs/output/SRR33542395/hisat/indx/GRCm39_index
```
-  hisat2-build - 	The HISAT2 utility that builds the index files from a reference genome (FASTA format) 
-  /home/vzscyborb... - path to the reference genome 
-  grcm_index - This is the base name for the output index files. HISAT2 will generate 8 files named grcm_index.1.ht2, grcm_index.2.ht2, ..., grcm_index.8.ht2.

### 2. **Alignment**

For the alignment require multiple files for 2 input method one for genome index and another is cleaned fastq file generated from fastp

```bash
hisat2 -p 1 \
       -x /home/vzscyborg/ngs/output/SRR33542395/hisat/indx/grcm_index \
       -U /home/vzscyborg/ngs/output/SRR33542395/fp/SRR33542395_clean.fastq \
       -S /home/vzscyborg/ngs/output/SRR33542395/hisat/align/SRR33542395_hisat.sam
```
-  hisat2 : HISAT2 is a fast and sensitive alignment tool for mapping sequencing reads (usually RNA-seq) to a reference genome. 
- -x : Specifies the basename of the HISAT2 index for the reference genome to which you want to align your reads. This is the index you built earlier from the mouse genome (GRCm39). 
-  -p 1 : Use 1 CPU thread for the alignment process. This limits HISAT2 to run on a single processor core 
-  -U : Specifies the input file containing unpaired single-end reads in FASTQ format. This is the cleaned sequencing reads you want to align. 
-  -S : Specifies the output file where the aligner will write the results in SAM format (Sequence Alignment/Map). This file contains detailed information on how each read aligns to the reference genome

---

## ▶️ Run Samtools

### 1. **Convert SAM → BAM**
```bash
samtools view -@ 3 -s 0.5 -b /home/vzscyborg/ngs/output/SRR33542395/hisat/align/SRR33542395_hisat.sam > /home/vzscyborg/ngs/output/SRR33542395/smtl/SRR33542395.bam
```
- samtools view : Converts between different formats (SAM ⇄ BAM), filters, and subsamples reads 
- -@ : 	Use 3 CPU threads to speed up the conversion. 
- -s 0.5 : 	Subsample 50% of reads randomly. The 0.5 means each read has a 50% chance of being retained. This is useful for downsampling large datasets  ( OPTIONAL )
- -b Output in BAM format (binary, compressed form of SAM). 
- 	Redirects output to a file.


### 2. **Sort BAM**
```bash
samtools sort -@ 3 -o /home/vzscyborg/ngs/output/SRR33542395/smtl/SRR35542395_sorted.bam /home/vzscyborg/ngs/output/SRR33542395/smtl/SRR33542395.bam
```

-  samtools sort : Convert Bam file into sorted bam file.
-  -@ : Use 3 CPU threads to speed up the conversion.
-  -o : specify the output file followed by the output file path along with the filename with .bam file extension
-  in last the directiory to the input bam file

### 3. **Index Sorted BAM**
```bash
samtools index -@ 3 -o /home/vzscyborg/ngs/output/SRR33542395/smtl/SRR33542395.bai /home/vzscyborg/ngs/output/SRR33542395/smtl/SRR33542395_sorted.bam
```

-  samtools index : Convert Sorted  Bam file into indexed bam file.
-  -@ : Use 3 CPU threads to speed up the conversion.
-  -o : specify the output file followed by the output file path along with the filename with .bam file extension
-  in last the directiory to the input sorted bam file

### **Chain all steps:**
```bash
samtools view -@ 1 -bS /home/vzscyborg/ngs/output/SRR33542395/hisat/align/SRR33542395_hisat.sam | \
samtools sort -@ 1 -o /home/vzscyborg/ngs/output/SRR33542395/smtl/SRR33542395_sorted.bam
samtools index -@ 1 /home/vzscyborg/ngs/output/SRR33542395/smtl/SRR33542395_sorted.bam
```

---

## ▶️ Run FeatureCounts

- For **gene/exon count matrix**, download reference annotation (e.g., from GENCODE) matching your genome build.
- Use `featureCounts` with sorted BAM and annotation files.

```
featureCounts -T 3 -t exon -g gene_id -a /home/vzscyborg/ngs/datasets/gtf/gencode.vM36.primary_assembly.annotation.gtf -o /home/vzscyborg/ngs/output/SRR33542395/fc/counts.txt /home/vzscyborg/ngs/output/SRR33542395/smtl/SRR33542395_sorted.bam
```

- featureCounts: program from the Subread package used to count reads mapped to genomic features (like genes or exons). 
- -T 3: 	Use 3 CPU threads for faster performance.  
- -t exon: 	Count exons as the feature type. Only lines with feature type == exon in the GTF file will be used. 
- -g gene_id: 	Group exons by the gene_id attribute from the GTF file to produce gene-level counts. 
- -a:  Path to the GTF annotation file, used to define the gene and exon features. 
- -o: 	Output file name where the read counts per gene will be saved.
- 	Input sorted BAM file (your aligned reads) that will be counted. 
---

<p align="center">
  <img src="https://capsule-render.vercel.app/api?type=waving&color=0:e96443,100:904e95&height=80&section=footer&animation=twinkling" />
</p>

<p align="center">
  <b>For more educational guides,contact https://www.linkedin.com/in/vidit-zainith-196960319?utm_source=share&utm_campaign=share_via&utm_content=profile&utm_medium=android_app.</b>
</p>
