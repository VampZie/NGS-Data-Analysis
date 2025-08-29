 # 🔗 Step 03: Read Alignment with BWA-MEM

## 🎯 Why Do We Do It?
Alignment maps sequencing reads to the reference genome to produce **BAM files** for downstream variant calling. `bwa mem` is widely used for short-read DNA alignment due to its speed and accuracy. The output BAM is compressed and indexed for efficient processing.

Use BWA-MEM when:
1. Working with Illumina short paired-end reads  
2. High-quality reference genome is available  
3. Preparing input for variant calling pipelines  

## 🛠️ Command

```bash
bwa mem -M -t "$th" "$ref" \
    "${d2}${sample}_1.fq" "${d2}${sample}_2.fq" | \
    samtools view -b -o "${d1}${sample}.bam" -
```

###🧾 Parameters Explained

```-M``` : Marks shorter split hits as secondary (compatible with Picard/GATK) 
```-t``` : Number of threads ($th specifies CPU cores) 
```$ref``` : Reference genome fasta file 
```${d2}${sample}_1.fq``` : Read 1 FASTQ file 
```${d2}${sample}_2.fq``` : Read 2 FASTQ file 
```samtools view -b``` : Converts SAM output to BAM format 
```-o``` : Output BAM file destination ${d1}${sample}.bam
