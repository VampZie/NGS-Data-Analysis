
# 🏷️ Step 04: Add or Replace Read Groups

## 🎯 Why Do We Do It?
Read groups are essential metadata tags required by **GATK** and other variant calling tools. They define the origin of reads (sample, platform, library, and sequencing run). Adding or replacing read groups ensures compatibility with downstream steps like BQSR and variant calling.

Use this step when:
1. BAM file has missing or incorrect read group information  
2. Preparing data for GATK pipelines  
3. Handling multiple samples, libraries, or sequencing runs  

## 🛠️ Command

```bash
for i in 34 35
do
  gatk AddOrReplaceReadGroups \
      -I  "${d1}${sample}.bam" \
      -O "${d1}${sample}.rg.bam" \
      -RGID  "${sample}" \
      -RGLB "lib1" \
      -RGPL "illumina" \
      -RGPU "${sample}.unit1" \
      -RGSM "${sample}"
done
```

### 🧾 Parameters Explained  
```-I``` : Input BAM file 

```-O``` : Output BAM file with updated read groups 

```-RGLB``` : Read Group Library (e.g., lib1) 

```-RGPL``` : Platform used (e.g., illumina) 

```-RGPU``` : Platform Unit (flowcell + lane info, here ${sample}.unit1) 

```-RGSM``` : Sample name (links reads to biological sample)
