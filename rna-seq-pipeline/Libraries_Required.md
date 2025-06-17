# 🧰 Required Tools & Libraries for RNA-Seq Analysis

This document outlines all the necessary software tools and libraries used in the RNA-Seq pipeline, along with installation instructions and usage notes. This setup is designed for Linux-based environments (Ubuntu/Debian).

---

## 📦 System-Level Dependencies

### `gunzip`
- **Purpose**: Decompress `.gz` compressed FASTQ files.
- **Install**:
  ```bash
  sudo apt install gzip
  ```
### build-essential, zlib1g-dev
Purpose: Compiler and development libraries required for building software like Subread.

Install: ```sudo apt install build-essential zlib1g-dev ```

### 🔄 Data Retrieval Tools
sra-tools
Purpose: Download sequencing data directly from NCBI using prefetch and fasterq-dump.

Install: ```sudo apt install sra-toolkit```

### 🧪 Quality Control Tools
FastQC
Purpose: Perform initial quality checks on raw FASTQ files.

Install: ```sudo apt install fastqc```

### ✂️ Read Trimming Tools
fastp
Purpose: Trim adapters, filter reads based on quality. Faster and more integrated than traditional trimmers.

Install:```sudo apt install fastp```

### 🎯 Alignment Tools
hisat2
Purpose: Fast and memory-efficient splice-aware aligner suitable for large genomes.

Install:```sudo apt install hisat2```

### 🧬 SAM/BAM File Processing
samtools
Purpose: View, convert, sort, and index alignment files (SAM/BAM/CRAM).

Install: ```sudo apt install samtools```

###  🧮 Feature Quantification
featureCounts (via Subread package)
Purpose: Assign mapped reads to genomic features such as genes or exons.

Install:

Download from SourceForge:
https://sourceforge.net/projects/subread/files/

Extract and build:
```
tar -xvzf subread-<version>.tar.gz
cd subread-<version>/src
make -f Makefile.Linux
```
(Optional) Add to PATH:
```
echo 'export PATH=$PATH:/path/to/subread/bin' >> ~/.bashrc
source ~/.bashrc
```
### 📊 Visualization & Summary Reports
MultiQC (optional but highly recommended)
Purpose: Aggregates output from FastQC, featureCounts, etc., into a single summary report.

Install:```pip install multiqc```


### 🧪 Differential Expression Analysis (R & Bioconductor)

R Base Installation
Install: ```sudo apt install r-base```

Bioconductor Packages in R
Launch R console: ```r```

Inside R:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("AnnotationDbi")
```
###✅ Check Installations
```
which fastqc
which hisat2
which samtools
which featureCounts
```

---
##⚠️ Notes
Some packages may require additional permissions or dependencies based on your system.

Keep your ~/.bashrc or ~/.zshrc updated with $PATH exports for ease of use.

Use Conda for a more portable and reproducible environment setup if sharing this pipeline across systems.
---
