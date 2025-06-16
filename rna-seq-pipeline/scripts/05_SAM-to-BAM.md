# 🔄 Step 05: Convert SAM to BAM Using Samtools

## 🎯 Why Do We Do It?

SAM (Sequence Alignment/Map) files are human-readable but **large and inefficient** for storage and processing. BAM (Binary Alignment/Map) files are the **compressed, binary equivalent**—ideal for downstream processing, sorting, and indexing. This step converts `.sam` to `.bam` format using `samtools`.
---



## 🛠️ Command

```
samtools view -@ 4  -bS /home/vzscyborg/ngs/mouse/hisat/31.sam > /home/vzscyborg/ngs/mouse/smtl/31.bam
```
