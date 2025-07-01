# 🧪 Step 03: Reference Genome Indexing with BWA

## 🎯 Why do we do it?

Before aligning sequencing reads to a reference genome, we need to **index the reference genome** using `bwa`. Indexing creates data structures that enable efficient and accurate alignment of reads. Without indexing, the aligner cannot map reads to the reference genome.

---

## 🛠️ Command

**Run reference genome indexing:**

```bash
bwa index ~/dsa/ref/GRCh38.primary_assembly.genome.fa
```
