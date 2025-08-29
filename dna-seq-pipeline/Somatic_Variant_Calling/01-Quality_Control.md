# 🧪 Step 01: Quality Control with FastQC

## 🎯 why do we do it?
Assess the quality of raw sequencing reads using **FastQC**, a tool that provides comprehensive metrics such as per-base quality, GC content, sequence duplication, and adapter contamination. This step ensures that only high-quality data proceed to downstream analysis.

---

## 🛠️ Command

**Run quality check:**
```bash
for i in 34 35
do
  fastqc -t 12 SRR76658${i}_1.fastq.gz SRR76658{i}_2.fastq.gz -O ~/fastq_qc_reports
done
```

- fastqc - command calling the quality check function
- -t refers to the threads using depanding upon the per core threads
- two reads input files for processing
- -O refers the path where the report of the quality check will save followed by the proper path
- SRRxxxxxxxxx.fastq.gz - files on which quality check will run
