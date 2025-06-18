# 🧬 RNA-Seq Pipeline for Mouse Transcriptomics

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

This repository provides a **modular, shell-and-R-based RNA-seq analysis workflow** for bulk transcriptome analysis in **mouse (Mus musculus)**. The pipeline is structured for clarity, educational reproducibility, and low-barrier adoption by students or researchers.

---

## 📁 Project Structure  

```
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
```


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


#### Sample Configuration :
[in `config/samples.csv`](https://github.com/VampZie/NGS-Data-Analysis/blob/b3a36de5e152222c6d54e06f289688cdacab4b15/rna-seq-pipeline/config/sample.csv)

---

## 🧬 Reference Files
### 🔸 Reference Genome (FASTA)
GRCm39.primary_assembly.genome.fa

Source: GENCODE Mouse

### 🔸 Annotation (GFF3)
gencode.vM37.primary_assembly.basic.annotation.gff3

Used in featureCounts for read quantification.

[config/paths.conf](https://github.com/VampZie/NGS-Data-Analysis/blob/7c5f0776c35780074f6a0f699f85ba4311f0289b/rna-seq-pipeline/config/paths.conf)



---
## ⚙️ Software & Environment
See [```/Libraries_Required.md```](https://github.com/VampZie/NGS-Data-Analysis/blob/cb604c83a31ba5f526db246dcb9e25d67d319523/rna-seq-pipeline/Libraries_Required.md) for installation details.

---

## 📊 Analysis Workflow
Each step is documented in [/script](https://github.com/VampZie/NGS-Data-Analysis/tree/a9ab6658ce3ef498704345eb79ba0f903220b381/rna-seq-pipeline/scripts) with detailed explanations and commands.

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
This project is licensed under the MIT License – see the [LICENSE](https://github.com/VampZie/NGS-Data-Analysis/blob/f9385362892099133b6d5d70bb84157fb688183b/LICENSE) file for details.

---
## 🧑 Author
Maintained by VampZie
Note: This is a personal, not institutional, GitHub identity.
---

## 💬 Feedback

Found a bug, typo, or want to suggest an improvement?  
Please [open an issue](https://github.com/VampZie/rna-seq-pipeline/issues) or submit a pull request.

📫 You can also reach out to me directly via [GitHub Discussions](https://github.com/VampZie/rna-seq-pipeline/discussions) or comments on this repository.
