# 🧬 RNA-Seq Pipeline for Mouse Transcriptomics

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

This repository provides a **modular, shell-and-R-based RNA-seq analysis workflow** for bulk transcriptome analysis in **mouse (Mus musculus)**. The pipeline is structured for clarity, educational reproducibility, and low-barrier adoption by students or researchers.

---

## 📁 Project Structure

rna-seq-pipeline/
├── config/                   # Configuration files
│   ├── paths.conf
│   ├── samples.csv
│   └── README.md
│
├── data/                     # Raw and processed data (NOT uploaded here)
│   └── README.md             # Instructions to obtain raw data
│
├── scripts/                  # Shell and R scripts for each analysis step
│   ├── 01_Quality-Control.md
│   ├── 02_Trimming.md
│   ├── 03_Genome-Indexing.md
│   ├── 04_Alignment.md
│   ├── 05_SAM-to-BAM.md
│   ├── 06_Sorting.md
│   ├── 07_BAM-Indexing.md
│   ├── 08_featureCount.md
│   └── 09_differential_gene_expression.md
│
├── bash_commands.md          # Optional: useful one-liners
├── Libraries_Required.md     # Required tools and installation
├── LICENSE                   # License file (MIT)
└── README.md                 # Main project documentation

---

## 🧪 Input Data

### 🔹 Experimental Design

- **HFD FMT (High-Fat Diet)**: 2 samples  
- **ND FMT (Normal Diet)**: 2 samples  
- Format: paired-end, mouse fecal microbiota transfer samples.

### 🔹 Data Sources

- Data downloaded from NCBI SRA using `prefetch` and `fasterq-dump`.
- Files are structured like:
data/raw/HFD1_1.fastq.gz
data/raw/HFD1_2.fastq.gz


#### Sample Configuration (in `config/samples.csv`):

```csv
sample,condition,fastq1,fastq2
HFD1,HFD,data/raw/HFD1_1.fastq.gz,data/raw/HFD1_2.fastq.gz
HFD2,HFD,data/raw/HFD2_1.fastq.gz,data/raw/HFD2_2.fastq.gz
ND1,ND,data/raw/ND1_1.fastq.gz,data/raw/ND1_2.fastq.gz
ND2,ND,data/raw/ND2_1.fastq.gz,data/raw/ND2_2.fastq.gz
```
---

## 🧬 Reference Files
### 🔸 Reference Genome (FASTA)
GRCm39.primary_assembly.genome.fa

Source: GENCODE Mouse

### 🔸 Annotation (GFF3)
gencode.vM37.primary_assembly.basic.annotation.gff3

Used in featureCounts for read quantification.

config/paths.conf example:
```
GENOME=/home/vzscyborg/ngs/mouse/GRCm39.primary_assembly.genome.fa
GTF=/home/vzscyborg/ngs/mouse/gencode.vM37.primary_assembly.basic.annotation.gff3
THREADS=8
```


---
## ⚙️ Software & Environment
See ```/Libraries_Required.md``` for installation details.

---

## 📊 Analysis Workflow
Each step is documented in /docs/ with detailed explanations and commands.

Quality Control – 01_Quality-Control.md

Trimming – 02_Trimming.md

Indexing Genome – 03_Genome-Indexing.md

Alignment – 04_Alignment.md

SAM to BAM – 05_SAM-to-BAM.md

Sorting – 06_Sorting.md

BAM Indexing – 07_BAM-Indexing.md

Quantification – 08_featureCount.md

Differential Expression – 09_differential_gene_expression.md

---

## 📌 Usage Notes
🧷 This repository is educational and modular — adapt to your needs.

🧼 Do not include raw FASTQ or BAM files in your repo. Use .gitignore.

📜 Scripts are meant to run step-by-step with clarity, not automation (for now).

---

## 📄 License
This project is licensed under the MIT License – see the LICENSE file for details.

---
## 🧑 Author
Maintained by VampZie
Note: This is a personal, not institutional, GitHub identity.
---

## 💬 Feedback

Found a bug, typo, or want to suggest an improvement?  
Please [open an issue](https://github.com/VampZie/rna-seq-pipeline/issues) or submit a pull request.

📫 You can also reach out to me directly via [GitHub Discussions](https://github.com/VampZie/rna-seq-pipeline/discussions) or comments on this repository.
