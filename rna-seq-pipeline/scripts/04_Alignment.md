khud h# 🧬 Step 04: Alignment Using HISAT2

## 🎯 Why Do We Do It?
Alignment is a **critical step** in RNA-seq and DNA-seq workflows where sequencing reads are mapped to a reference genome. This allows us to identify where each read originated from, enabling downstream analysis like gene expression quantification and variant discovery.

We use **HISAT2** for this task due to its high speed, low memory footprint, and support for spliced alignment, which is essential for RNA-seq.
---



## 🛠️ Command 

for file having single read
```bash 
hisat2 -p 1 \
       -x /home/vzscyborg/ngs/output/SRR33542395/hisat/indx/grcm_index \
       -U /home/vzscyborg/ngs/output/SRR33542395/fp/SRR33542395_clean.fastq \
       -S /home/vzscyborg/ngs/output/SRR33542395/hisat/align/SRR33542395_hisat.sam
```

🧾 Parameters Explained
Parameter |	Description
---
hisat2	Command to launch the aligner
-p 1	Use 1 CPU thread for alignment
-x	Path to the HISAT2 genome index (built in Step 3)
-U	Input file for single-end reads (FASTQ format)
-S	Output file for alignments in SAM format

for Data having two reads
```bash 
hisat2 -p 4 -x /home/vzscyborg/ngs/datasets/grcm39/m39_index  fastq -1 /home/vzscyborg/ngs/mouse/fp/31r1c.fastq -2 /home/vzscyborg/ngs/mouse/fp/31r2c.fastq -S /home/vzscyborg/ngs/mouse/hisat/31.sam
```

-  hisat2 : HISAT2 is a fast and sensitive alignment tool for mapping sequencing reads (usually RNA-seq) to a reference genome. 
- -x : Specifies the basename of the HISAT2 index for the reference genome to which you want to align your reads. This is the index you built earlier from the mouse genome (GRCm39). 
-  -p 1 : Use 1 CPU thread for the alignment process. This limits HISAT2 to run on a single processor core 
-  -U : Specifies the input file containing unpaired single-end reads in FASTQ format. This is the cleaned sequencing reads you want to align. 
-  -S : Specifies the output file where the aligner will write the results in SAM format (Sequence Alignment/Map). This file contains detailed information on how each read aligns to the reference genome
