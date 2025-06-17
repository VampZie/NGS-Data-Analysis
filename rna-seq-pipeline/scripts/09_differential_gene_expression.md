# 📈 Step 09: Differential Gene Expression Analysis (DGE)

## 🎯 Why Do We Do It?
Differential Gene Expression analysis identifies genes whose expression levels significantly differ between experimental conditions (e.g., treated vs. control). This is essential for understanding biological responses and regulatory mechanisms.

---

## 📦 Tools Required
This analysis is typically performed in **R** using one of the following packages:
- `DESeq2` – Recommended for count data from RNA-seq
- `edgeR` – Ideal for datasets with biological replicates
- `limma-voom` – Best for large-scale datasets

---

## 🛠️ Workflow Overview

### 1. **Prepare Input**
- Use the `counts.txt` file generated from `featureCounts` (Step 08).
- Prepare a **sample metadata** file (CSV or TSV) with columns like: `sample`, `condition`, `replicate`, `file_path`. In my case i used 4 samples 2 HFD FMT and 2 ND FMT , so i manually placed it in the code of ```R```.
Example:
```csv
sample,condition,replicate,file
Treated_1,Treated,1,/path/to/Treated_1.bam
Treated_2,Treated,2,/path/to/Treated_2.bam
Control_1,Control,1,/path/to/Control_1.bam
Control_2,Control,2,/path/to/Control_2.bam
```
---

## R Code

## 📘 R Code Disclaimer

This R script includes a range of functions, parameters, and configurations that may not all be relevant to your specific use case.

Please consider the following before using:

- 🔧 **Modify accordingly**: Not all features may apply to your dataset or computational environment. Adjust the code to suit your needs.
- ✅ **Use responsibly**: Ensure you're executing the code in an appropriate and authorized environment.
- 🛡️ **Respect resource limits**: Avoid overloading shared or institutional computational resources. Unauthorized or excessive usage may lead to restrictions or loss of access.

> ⚠️ **IMPORTANT**: All code provided is intended for educational and research purposes. Use it ethically, and always cite relevant tools or libraries when publishing results.



```
library(DESeq2)

counts <- read.table("~/ngs/mouse/counts.txt", header=TRUE, row.names=1, check.names=FALSE)
# Drop annotation columns if present
counts <- counts[ , 6:ncol(counts) ]  # adjust if there are fewer columns before samples

# Manually define sample conditions
condition <- factor(c("HFD", "HFD", "ND", "ND"))  # order must match column order
coldata <- data.frame(row.names = colnames(counts), condition)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)

# Run DESeq2
dds <- DESeq(dds)

# Extract and save results
res <- results(dds)
resOrdered <- res[order(res$padj), ]
write.csv(as.data.frame(resOrdered), file="deseq2_results.csv")



#Filter Significant Genes
# Adjust p-value threshold (commonly 0.05 or 0.01)
sig_genes <- res[which(res$padj < 0.05), ]

# Optional: Add log2 fold change threshold
sig_genes <- sig_genes[abs(sig_genes$log2FoldChange) > 1, ]

# Save
write.csv(as.data.frame(sig_genes), "sig_genes.csv")
# Adjust p-value threshold (commonly 0.05 or 0.01)
sig_genes <- res[which(res$padj < 0.05), ]

# Optional: Add log2 fold change threshold
sig_genes <- sig_genes[abs(sig_genes$log2FoldChange) > 1, ]

# Save
write.csv(as.data.frame(sig_genes), "sig_genes.csv")


#MA Plot (for gene-wise expression shift)
plotMA(res, ylim=c(-5,5))


#Volcano Plot
library(ggplot2)

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df$significant <- res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1

ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
  geom_point(alpha=0.4) +
  scale_color_manual(values=c("grey", "red")) +
  theme_minimal() +
  labs(title="Volcano Plot", x="Log2 Fold Change", y="-log10 adjusted p-value")



# Gene Annotation
BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
library(AnnotationDbi)

gene_symbols <- mapIds(org.Mm.eg.db, 
                       keys = rownames(sig_genes),
                       column = "SYMBOL",                 
                       keytype = "ENSEMBL",
                       multiVals = "first")

sig_genes$symbol <- gene_symbols
write.csv(as.data.frame(sig_genes), "sig_genes_annotated.csv")


# Load libraries
library(DESeq2)
library(org.Mm.eg.db)
library(AnnotationDbi)

# Read featureCounts outputs
fc1 <- read.table("HFD1.txt", header=TRUE, sep="\t", comment.char="#", row.names=1)
fc2 <- read.table("HFD2.txt", header=TRUE, sep="\t", comment.char="#", row.names=1)
fc3 <- read.table("ND1.txt",  header=TRUE, sep="\t", comment.char="#", row.names=1)
fc4 <- read.table("ND2.txt",  header=TRUE, sep="\t", comment.char="#", row.names=1)

# Extract only the counts (7th column onwards)
count_data <- data.frame(
  HFD1 = fc1[, 7],
  HFD2 = fc2[, 7],
  ND1  = fc3[, 7],
  ND2  = fc4[, 7]
)

# Clean rownames to remove Ensembl version suffix
rownames(count_data) <- sub("\\..*", "", rownames(count_data))

# Build sample info
col_data <- data.frame(
  row.names = colnames(count_data),
  condition = c("HFD", "HFD", "ND", "ND")
)

# Build DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = col_data,
                              design = ~ condition)

# Run DESeq2
dds <- DESeq(dds)
res <- results(dds)

# Clean Ensembl IDs again for mapping
clean_ids <- sub("\\..*", "", rownames(res))

# Map to gene symbols
symbols <- mapIds(org.Mm.eg.db,
                  keys = clean_ids,
                  column = "SYMBOL",
                  keytype = "ENSEMBL",
                  multiVals = "first")

# Add gene symbols to results
res$symbol <- symbols

# Sort by adjusted p-value
res_ordered <- res[order(res$padj), ]

# Save to file
write.csv(as.data.frame(res_ordered), file = "deseq2_results.csv")

# Print summary
summary(res)
```
---
---
## 🧠 Interpretation
Look for:

Genes with adjusted ```p-value``` < ```0.05```

```Log2 Fold``` Change thresholds (e.g., >1 or <-1) for biological significance

Pathway enrichment can follow using tools like ```clusterProfiler```, ```GSEA```, or ```Enrichr```
---
---
