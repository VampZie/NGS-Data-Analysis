# ✂️ Step 02: Run Cutadapt if FastQC Flags Issues

## 🎯 Why Do We Do It?
If the **FastQC report reveals quality issues**, `cutadapt` is used for **adapter trimming**, **quality filtering**, and **length control** in paired-end sequencing reads. It allows explicit adapter specification, removes low-quality bases, and filters short reads.

Use Cutadapt **when the following conditions arise** in FastQC reports:
1. Adapter contamination  
2. Low-quality bases at 5'/3' ends  
3. Biased base content at read ends  
4. Overrepresented adapter dimers or primers  
5. Reads shorter than a usable threshold  

## 🛠️ Command

```bash

for i in 34 35
do 
 cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
          -A AGATCGGAAGAGGCGTCGTGTAGGGAAAGAGTT \
          -o ${d2}${i}_1.fastq \
          -p ${{d2}${i}_2.fastq \
          ${d0}34_1.fastq \
          ${d0}34_2.fastq \
          --minimum-length 15 \
          --cores=12
```

### 🧾 Parameters Explained 
```-a``` : Adapter sequence for Read 1 ( NOTE: adapter sequence may varies according to your datasets please refer to the Fastqc report ) 
```-A``` : Adapter sequence for Read 2 ( NOTE: adapter sequence may varies according to your datasets please refer to the Fastqc report ) 
```-o``` : Output file for trimmed Read 1 
```-p``` : Output file for trimmed Read 2 
```--minimum-length``` : Discard reads shorter than threshold 
```--cores``` : Number of threads to use 
