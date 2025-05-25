<!-- HEADER BANNER -->
<p align="center">
  <img src="https://capsule-render.vercel.app/api?type=waving&color=0:e96443,100:904e95&height=150&section=header&text=RNA-Seq%20Preparation%20&%20FASTQ%20Download%20Guide&fontSize=32&fontColor=fff&animation=twinkling" alt="Header Banner"/>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/Educational%20Resource-Bioinformatics-blueviolet?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/Level-Beginner%20to%20Intermediate-00bfff?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/Workflow-RNA--Seq%20Pipeline-ff69b4?style=for-the-badge"/>
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

For NCBI, use the **SRA Toolkit**.

- Check install:  
  ```bash
  vdb-dump --version
  ```

- Download using `fasterq-dump`:
  ```bash
  fasterq-dump SRR33542395 --split-files -O ~/Downloads --temp ~/Downloads/tmp --progress
  ```
  - `--split-files` : Downloads both reads (if paired-end)
  - `-O` : Output directory
  - `--temp` : Directory for temp files (useful for large downloads)
  - `--progress` : Shows download progress

---

# 🧬 STARTING RNA-SEQ ANALYSIS

---

## ▶️ Run FastQC on the Downloaded Files

**Run quality check:**
```bash
fastqc SRR1551114_1.fastq.gz SRR1551114_2.fastq.gz -o ~/fastq_qc_reports
# or
fastqc -t 4 -o /home/vzscyborg/rnaseq/fastqc_rep /home/vzcyborg/datasets/rnads/SRR1551114_1.fastq.gz
```
- `-t 4` : Number of threads
- `-o` : Output directory
- `SRRxxxxxxx.fastq.gz` : Input files

---

## ▶️ Run Fastp, if FastQC Report Fails

Run **Fastp** for:
  1. Adapter contamination
  2. Low-quality tiles
  3. Biased base content
  4. Overrepresented sequences

```bash
fastp -i /home/vzscyborg/ngs/datasets/SRR33542395 \
      -o /home/vzscyborg/ngs/output/SRR33542395/fp/SRR33542395_clean.fastq \
      --detect_adapter_for_pe \
      --html /home/vzscyborg/ngs/output/SRR33542395/fp/SRR33542395_clean.html \
      --json /home/vzscyborg/ngs/output/SRR33542395/fp/SRR33542395_clean.json \
      --thread 1
```
- `-i` : Input file (for paired-end, use `-I` for Read 2)
- `-o` : Output file (for paired-end, use `-O` for Read 2)
- `--html` : Output report (HTML)
- `--thread` : Number of threads

---

## ▶️ Run HISAT2

### 1. **Genome Indexing**
Download your reference genome (e.g., GRCm39 from NCBI) in FASTA format.

```bash
hisat2-build /home/vzscyborg/ngs/datasets/gtf/GRCm39.primary_assembly.genome.fa \
             /home/vzscyborg/ngs/output/SRR33542395/hisat/indx/GRCm39_index
```
- `hisat2-build` : Create index from FASTA
- Output: GRCm39_index.*.ht2

### 2. **Alignment**

```bash
hisat2 -p 1 \
       -x /home/vzscyborg/ngs/output/SRR33542395/hisat/indx/grcm_index \
       -U /home/vzscyborg/ngs/output/SRR33542395/fp/SRR33542395_clean.fastq \
       -S /home/vzscyborg/ngs/output/SRR33542395/hisat/align/SRR33542395_hisat.sam
```
- `-x` : Path to genome index
- `-U` : Input clean reads
- `-S` : Output SAM file

---

## ▶️ Run Samtools

### 1. **Convert SAM → BAM**
```bash
samtools view -@ 3 -s 0.5 -b /home/vzscyborg/ngs/output/SRR33542395/hisat/align/SRR33542395_hisat.sam > /home/vzscyborg/ngs/output/SRR33542395/smtl/SRR33542395.bam
```
- `-@ 3` : Threads
- `-s 0.5` : Subsample 50% (optional)
- `-b` : Output BAM

### 2. **Sort BAM**
```bash
samtools sort -@ 3 -o /home/vzscyborg/ngs/output/SRR33542395/smtl/SRR35542395_sorted.bam /home/vzscyborg/ngs/output/SRR33542395/smtl/SRR33542395.bam
```

### 3. **Index Sorted BAM**
```bash
samtools index -@ 3 -o /home/vzscyborg/ngs/output/SRR33542395/smtl/SRR33542395.bai /home/vzscyborg/ngs/output/SRR33542395/smtl/SRR33542395_sorted.bam
```

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

---

<p align="center">
  <img src="https://capsule-render.vercel.app/api?type=waving&color=0:e96443,100:904e95&height=80&section=footer&animation=twinkling" />
</p>

<p align="center">
  <b>For more educational guides, visit <a href="https://www.ebi.ac.uk/training/online/courses">EBI Training</a> or <a href="https://rnabio.org">RNA-Seq Tutorials</a>.</b>
</p>
