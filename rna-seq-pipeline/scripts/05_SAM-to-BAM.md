# 🔄 Step 05: Convert SAM to BAM Using Samtools

## 🎯 Why Do We Do It? 


SAM (Sequence Alignment/Map) files are human-readable but **large and inefficient** for storage and processing. BAM (Binary Alignment/Map) files are the **compressed, binary equivalent**—ideal for downstream processing, sorting, and indexing. This step converts `.sam` to `.bam` format using `samtools`.


---



## 🛠️ Command

```
samtools view -@ 4  -bS /home/vzscyborg/ngs/mouse/hisat/31.sam > /home/vzscyborg/ngs/mouse/smtl/31.bam
```

### 🧾 Parameters

- samtools view : Converts between different formats (SAM ⇄ BAM), filters, and subsamples reads 
- -@ : 	Use 3 CPU threads to speed up the conversion. 
- -s 0.5 : 	Subsample 50% of reads randomly. The 0.5 means each read has a 50% chance of being retained. This is useful for downsampling large datasets  ( OPTIONAL )
- -b Output in BAM format (binary, compressed form of SAM). 
- 	Redirects output to a file.


---

## ✅ Output 
This file will be used in the next steps for sorting, indexing, and feature counting.
