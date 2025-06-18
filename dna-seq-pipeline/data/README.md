# 🧬 DNA-seq Analysis Workflow – Matched Tumor-Normal from Human Prostate Cancer

This repository accompanies a comprehensive **DNA-seq analysis workflow** for human cancer, focusing on whole genome sequencing (WGS) of matched **tumor and normal** samples from **prostate adenocarcinoma**. It supports somatic variant calling, tumor evolution exploration, and potential driver mutation detection.

---

## 📁 Dataset Description

### 🔹 Biological Context

- **Organism:** *Homo sapiens* (Human)
- **Disease:** Prostate Adenocarcinoma
- **Study Focus:** Matched primary tumor and metastatic lymph node samples from a single patient, with matched normal tissue (blood).
- **Total Samples:** 5
  - `19651_Blood` – [Normal (Blood)](https://www.ncbi.nlm.nih.gov/biosample/14209982)
  - `19651_LP` – [Tumor (Left prostate biopsy)](https://www.ncbi.nlm.nih.gov/biosample/14209984)
  - `19651_RP` – [Tumor (Right prostate biopsy)](https://www.ncbi.nlm.nih.gov/biosample/14209986)
  - `19651_LLN` – [Metastatic lymph node (Left)](https://www.ncbi.nlm.nih.gov/biosample/14209983)
  - `19651_RLN` – [Metastatic lymph node (Right)](https://www.ncbi.nlm.nih.gov/biosample/14209985)

---

## 📥 Data Acquisition

### 🔹 1. Raw Sequencing Data

- **Source:** [NCBI SRA]([https://www.ncbi.nlm.nih.gov/sra](https://www.ncbi.nlm.nih.gov/biosample?LinkName=bioproject_biosample_all&from_uid=608841))
- **Tools Used:**
  - [`prefetch`](https://github.com/ncbi/sra-tools) – to download `.sra` files
  - [`fasterq-dump`](https://github.com/ncbi/sra-tools) – to convert `.sra` to `.fastq`

```bash
# Example commands
prefetch SRS6218633  # 19651_Blood
fasterq-dump SRS6218633
```
