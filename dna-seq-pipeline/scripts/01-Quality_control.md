# 🧪 Step 01: Quality Control with FastQC

## 🎯 why do we do it?
Assess the quality of raw sequencing reads using **FastQC**, a tool that provides comprehensive metrics such as per-base quality, GC content, sequence duplication, and adapter contamination. This step ensures that only high-quality data proceed to downstream analysis.

---

## 🛠️ Command

**Run quality check:**
```bash
# for Normal 
fastqc -t 4 SRR7665835_1.fastq SRR7665835_2.fastq -O ~/dsa/35/fqc --progress

# for Tumor bulk
fastqc -t 4 SRR7665834_1.fastq SRR7665834_2.fastq -O ~/dsa/34/fqc --progress

```
- fastqc - command calling the quality check function
- -t 4 - t refers to the threads using depanding upon the per core threads
- two reads input files for processing
- -O refers the path where the report of the quality check will save followed by the proper path
- SRRxxxxxxxxx.fastq.gz - files on which quality check will run
