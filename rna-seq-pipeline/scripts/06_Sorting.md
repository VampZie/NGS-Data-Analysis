# 📊 Step 06: Sort BAM File Using Samtools

## 🎯 Why Do We Do It? 

After converting a SAM file to BAM, the next critical step is to **sort the BAM file by genomic coordinates**. This is essential for downstream tools such as indexing, duplicate marking, and read counting, all of which require sorted input.

---

### 🔧 Tool: `samtools sort`

---

## 🛠️ Command 

```bash
samtools sort -@ 4 -o /home/vzscyborg/ngs/mouse/smtl/31sorted.bam /home/vzscyborg/ngs/mouse/smtl/31.bam
```
### 🧾 Parameters
-  samtools sort : Convert Bam file into sorted bam file.
-  -@ : Use 3 CPU threads to speed up the conversion.
-  -o : specify the output file followed by the output file path along with the filename with .bam file extension
-  in last the directiory to the input bam file


## ✅ Output

A sorted BAM file:
/home/vzscyborg/ngs/mouse/smtl/31sorted.bam

This sorted file is required for BAM indexing and gene-level feature counting in upcoming steps.
