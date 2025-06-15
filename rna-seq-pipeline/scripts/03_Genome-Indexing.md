# 🧬 Step 03: Genome Indexing Using HISAT2

## 🎯 Why Do We Do It?
Genome indexing is a **prerequisite for aligning sequencing reads** to a reference genome. It transforms the reference FASTA file into a format that enables **efficient and accurate alignment** by tools like HISAT2. Indexing significantly accelerates the alignment process while conserving memory during execution.

---



## 🔧 Tool: `hisat2-build`

The `hisat2-build` utility constructs a **hierarchical FM index** from a reference genome in FASTA format.

### Installation: 
``` sudo apt install hisat2 ```

---



## 📥 Reference Genome Used
 Species: Mus musculus (mouse) 

 Assembly: GRCm39 
 File Format: FASTA  
 
 File: GRCm39.primary_assembly.genome.fa ( [downloaded from GENCODE](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M37/GRCm39.primary_assembly.genome.fa.gz)  )
 

---



## 🛠️ Command

```bash
hisat2-build -p 2 /home/vzscyborg/ngs/datasets/gtf/GRCm39.primary_assembly.genome.fa \
             /home/vzscyborg/ngs/output/SRR33542395/hisat/indx/GRCm39_index
```
-  hisat2-build - 	The HISAT2 utility that builds the index files from a reference genome (FASTA format)
-  -p 2: number threads of the cpu using by this command 
-  /home/vzscyborb... - path to the reference genome 
-  grcm_index - This is the base name for the output index files. HISAT2 will generate 8 files named grcm_index.1.ht2, grcm_index.2.ht2, ..., grcm_index.8.ht2.

---



### 🔖 HISAT2 will generate eight index files:
GRCm39_index.1.ht2  
GRCm39_index.2.ht2  
...  
GRCm39_index.8.ht2  

---
---
---

This indexing step is critical and only needs to be performed once per reference genome. Once the index is built, it can be reused across multiple RNA-seq or DNA-seq alignment workflows using HISAT2.

---
---
---
