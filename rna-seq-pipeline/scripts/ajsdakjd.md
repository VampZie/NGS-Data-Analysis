# 🧪 Step 01: Quality Control with FastQC

## 🎯 why do we do it?
Assess the quality of raw sequencing reads using **FastQC**, a tool that provides comprehensive metrics such as per-base quality, GC content, sequence duplication, and adapter contamination. This step ensures that only high-quality data proceed to downstream analysis.

---

## 🛠️ Command

**Run quality check:**
```bash
fastqc SRR1551114_1.fastq.gz SRR1551114_2.fastq.gz -O ~/fastq_qc_reports --progress
```
or
```
fastqc -t 4 /home/vzscyborg/ngs/mouse/31r1.fastq  /home/vzscyborg/ngs/mouse/31r2.fastq -O /home/vzscyborg/ngs/mouse/fqc --progress
```
- fastqc - command calling the quality check function
- -t 4 - t refers to the threads using depanding upon the per core threads
- two reads input files for processing
- -O refers the path where the report of the quality check will save followed by the proper path
- SRRxxxxxxxxx.fastq.gz - files on which quality check will run
