----------
# WORK IS ONGOING
-----------

```
dna-seq-pipeline/
├── config/
│   ├── samples.csv              # tumor-normal metadata
│   ├── paths.conf               # reference genome, tool paths
│   └── README.md
│
├── data/                        # Instructions only, no uploads
│   └── README.md
│
├── scripts/
│   ├── 01_Download_FASTQ.md
│   ├── 02_Quality_Control.md
|   ├── 03_Trimming.md
|   ├── 04_Reference_genome_indexing.md
│   ├── 05_Alignment_BWA.md
│   ├── 06_Sort_Index_BAM.md
│   ├── 07_Mark_Duplicates.md
│   ├── 08_BaseRecalibration.md
│   ├── 09_Somatic_Variant_Calling.md   # GATK Mutect2
│   ├── 10_Filtering_Variants.md
│   ├── 11_Annotation.md                # VEP or ANNOVAR
│   ├── 12_CopyNumber_Analysis.md       # CNVkit or GATK CNV
│   └── 13_Clonality_Heterogeneity.md   # PyClone/SciClone
│
├── Libraries_Required.md
├── LICENSE
└── README.md
```
