# 🧬 DNA-seq Analysis Workflow – Matched Tumor-Normal from Her2+ Breast Cancer

This repository contains a comprehensive **DNA-seq analysis workflow** for human cancer, focusing on **whole-genome sequencing (WGS)** of **matched tumor and normal samples** from a Her2-positive breast cancer patient. The pipeline supports somatic variant calling, clonality analysis, and tumor heterogeneity assessment.

---

## 📁 Dataset Description

### 🔹 Biological Context

- **Organism:** *Homo sapiens* (Human)
- **Disease:** Her2-positive Breast Cancer
- **Study Focus:** Tumor cell clusters and matched normal tissue from a single patient
- **Total Samples:** 7  
  - [`SRR7665835`](https://www.ebi.ac.uk/ena/browser/view/SRR7665835) – **Normal** (blood/tissue matched)  
  - [`SRR7665834`](https://www.ebi.ac.uk/ena/browser/view/SRR7665834) – **Tumor Bulk**  
  - [`SRR7665833`](https://www.ebi.ac.uk/ena/browser/view/SRR7665833) – **Cell Cluster 1**  
  - [`SRR7665832`](https://www.ebi.ac.uk/ena/browser/view/SRR7665832) – **Cell Cluster 2**  
  - [`SRR7665831`](https://www.ebi.ac.uk/ena/browser/view/SRR7665831) – **Cell Cluster 3**  
  - [`SRR7665830`](https://www.ebi.ac.uk/ena/browser/view/SRR7665830) – **Cell Cluster 4**  
  - [`SRR7665829`](https://www.ebi.ac.uk/ena/browser/view/SRR7665829) – **Cell Cluster 5**  
---

## 📥 Data Acquisition

### 🔹 1. Raw Sequencing Data

- **Source:** [European Nucleotide Archieve]((https://www.ebi.ac.uk/ena/browser))  
- **Tools Required:**
  - wget
  - gunzip
#### 📦 Example Commands

```bash
# Download and convert Normal Sample
wget https://www.ebi.ac.uk/ena/browser/view/SRR7665835
gunzip SRR7665835.fastq.gz

# (Repeat similarly for other SRR IDs)
