# 🧬 RNA-seq Analysis Dataset Guide

This repository accompanies a basic-to-intermediate RNA-seq analysis workflow focused on mouse gene expression. It includes practical steps starting from raw read retrieval to differential gene expression (DGE) analysis using tools like `HISAT2`, `samtools`, `featureCounts`, and `DESeq2`.

---

## 📁 Dataset Description i used for this RNA-seq Analysis

### 🔹 Biological Context
- Organism: **Mus musculus** (Mouse)
- Experimental design: **FMT (Fecal Microbiota Transplantation)** under two conditions:
  - **HFD** – High-Fat Diet FMT (2 samples)
  - **ND** – Normal Diet FMT (2 samples)

---

## 📥 Data Acquisition

### 🔹 1. **Raw Sequencing Data**
- **Source:** [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra)
- **Tools used:**
  - [`prefetch`](https://github.com/ncbi/sra-tools) — for downloading `.sra` files
  - [`fasterq-dump`](https://github.com/ncbi/sra-tools) — for converting `.sra` to `.fastq`
  
```bash
# Download raw SRA files (replace with actual accession numbers)
prefetch SRRxxxxxxx

# Convert to FASTQ
fasterq-dump SRRxxxxxxx
```

Note: Ensure SRA Toolkit is configured with correct cache and output directories.

---

### 🔹 2. Reference Genome and Annotation

Both genome and annotation were downloaded from the **GENCODE** database for mouse:

#### 🔸 Genome FASTA  
- **File:** `GRCm39.primary_assembly.genome.fa`  
- **URL:** [https://www.gencodegenes.org/mouse/](https://www.gencodegenes.org/mouse/)

#### 🔸 Annotation GFF3  
- **File:** `gencode.vM37.primary_assembly.basic.annotation.gff3`  
- **Used for:** Read counting using `featureCounts`

---

### ⚠️ Notes

- The annotation format (`.gff3`) **must match** the genome build used during alignment.
- Always verify **chromosome naming consistency** (e.g., `chr1` vs `1`) between the **FASTA** and **GFF3/GTF** files.
- The `fasterq-dump` output will typically contain **paired-end files**:
  - `*_1.fastq`
  - `*_2.fastq`
 
---
