# ✂️ Step 02: Run Fastp if FastQC Flags Issues

## 🎯 Why Do We Do It?
If the **FastQC report reveals quality issues**, `fastp` is used for **adapter trimming**, **quality filtering**, and **base correction** in paired-end sequencing reads. It is a fast and versatile tool that can automatically detect adapters, remove low-quality reads, and generate informative quality reports.

Use Fastp **when the following conditions arise** in FastQC reports:
1. Adapter contamination  
2. Low-quality tiles or per-base quality drop  
3. Biased base content at 5'/3' ends  
4. Overrepresented sequences  
5. Poor overall read quality needing cleanup

## 🛠️ Command

```bash
fastp -i /home/vzscyborg/ngs/mouse/31r1.fastq  -I /home/vzscyborg/ngs/mouse/31r2.fastq -o /home/vzscyborg/ngs/mouse/fp/31r1c.fastq -O /home/vzscyborg/ngs/mouse/fp/31r2c.fastq --detect_adapter_for_pe   --html /home/vzscyborg/ngs/mouse/fp/31c.html pe   --json /home/vzscyborg/ngs/mouse/fp/31c.json --thread  4
```
### 🧾 Parameters Explained
- -i: for input file followed by file directory, if two read second read start with \-I (source directory) 
-  -o: for output file followed by file directory with extension, if two out reads second read start with /-O (destination directory) 
-  --html : for making the output report in html formate along with the file name 
-  --thread define the number of thread fastp should use

## ⚠️ Force Trimming (If Quality Issues Persist)
Even after Fastp, if reads fail quality checks or show insert size bias, force trimming and cleanup must be performed manually.

🔹 1. Filter Short Inserts
   Removes fragments shorter than 50 bp to eliminate noise from short, non-informative reads.
   ```
   fastp -i SRR33583544_1_clean.fastq -I SRR33583544_2_clean.fastq -o r1.fastq -O r2.fastq --length_required 50
   ```
   
🔹 2. Mark and Remove Duplicates
   Sort and deduplicate aligned reads to avoid bias in expression quantification.
   ```
   samtools sort -o sorted.bam input.bam
   samtools markdup sorted.bam dedup.bam
   ```
🔹 3. Check Insert Size Distribution
   Use Picard to confirm if the library has excessive short inserts or other anomalies.
   Confirm if the entire library is composed of short inserts:
   ```
   picard CollectInsertSizeMetrics I=dedup.bam O=insert_size_metrics.txt H=insert_size_histogram.pdf
   ```
   
