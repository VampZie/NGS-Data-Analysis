# NGS
All about Basic NGS
----------------------------

Data mining using RNA seq data 
==============================

1. Data Preprocessing
---------------------
Before mining, RNA-Seq data needs to be processed:

    + Quality Control (QC): Use tools like FastQC.

    + Trimming: Remove adapters/low-quality reads (e.g., Trimmomatic).

    + Alignment: Map reads to the genome (e.g., HISAT2, STAR).

    + Quantification: Count reads per gene/transcript (e.g., featureCounts, HTSeq).

    + Normalization: Use methods like TPM, RPKM, or DESeq2's normalization.


2. Feature Extraction
---------------------

Turn raw counts into a data matrix:
    + Rows = Genes
    + Columns = Samples
    + Values = Normalized expression levels
This matrix is your basis for downstream mining

3. Common Data Mining Tasks
---------------------------
Popular techniques used on RNA-Seq datasets:

	+ Differential Gene Expression (DGE)
    		Identify genes that are up- or down-regulated under different conditions.
    		Tools: DESeq2, edgeR, limma.

	+ Clustering
    		Group genes or samples based on expression profiles.
    		Techniques: K-means, Hierarchical clustering, DBSCAN.
    		Purpose: Discover co-expressed gene modules or sample subtypes.

	+ Classification
    		Predict sample labels (e.g., disease vs. control).
    		Methods: Random Forest, SVM, Neural Networks.
    		Requires labeled datasets.

	+ Dimensionality Reduction
    		Reduce complexity while retaining key information.
    		Methods: PCA, t-SNE, UMAP.
    		Used for visualization or noise reduction.

	+ Gene Co-expression Network Analysis
    		Construct networks showing how genes co-vary.
    		Tool: WGCNA (Weighted Gene Co-expression Network Analysis).
    		Helps find gene modules associated with traits.

	+ Functional Enrichment & Pathway Analysis
    		After finding differentially expressed genes or modules:
        	Perform GO enrichment, KEGG pathway analysis, etc.
        	Tools: clusterProfiler, gProfiler, Enrichr.

4. Advanced Techniques (for data mining research)
-------------------------------------------------
+ Deep Learning: Autoencoders, CNNs, RNNs for feature extraction or prediction.
+ Multi-Omics Integration: Combine RNA-Seq with methylation, proteomics, etc.
+ Text Mining: Link gene findings to literature using NLP.

5. Toolkits & Environment
    + R/Bioconductor: DESeq2, edgeR, limma, WGCNA.
    + Python: Scanpy, scikit-learn, pandas, TensorFlow.
    + Galaxy: For GUI-based workflows without coding.

=================================================================
-----------------------------------------------------------------
=================================================================

🔁 Typical RNA-Seq Mining Pipeline

    Get data (e.g., from GEO or SRA)

    QC, alignment, quantification

    Build expression matrix

    Normalize and filter

    Perform mining tasks (clustering, classification, etc.)

    Validate & interpret biologically
