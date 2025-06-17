# 🧷 Step 07: Index Sorted BAM File Using Samtools

## 🎯 Why Do We Do It?
Indexing a sorted BAM file enables **rapid random access** to alignment data. This is crucial for visualization tools (like IGV), variant calling, and read quantification, where only specific regions of the genome are queried.

---
## 🔧 Tool: `samtools index`
---

```bash
samtools index -@ 4 -o /home/vzscyborg/ngs/mouse/smtl/31.bai /home/vzscyborg/ngs/mouse/smtl/31sorted.bam
```

### 🧾 Parameters 
-  samtools index : Convert Sorted  Bam file into indexed bam file.
-  -@ : Use 3 CPU threads to speed up the conversion.
-  -o : specify the output file followed by the output file path along with the filename with .bam file extension
-  in last the directiory to the input sorted bam file

---
---
## ✅ Output
An indexed BAM file (.bai extension):
/home/vzscyborg/ngs/mouse/smtl/31.bai

This index file allows efficient data retrieval during visualization or transcript quantification.
