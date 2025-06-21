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
│   ├── 03_Alignment_BWA.md
│   ├── 04_Sort_Index_BAM.md
│   ├── 05_Mark_Duplicates.md
│   ├── 06_BaseRecalibration.md
│   ├── 07_Somatic_Variant_Calling.md   # GATK Mutect2
│   ├── 08_Filtering_Variants.md
│   ├── 09_Annotation.md                # VEP or ANNOVAR
│   ├── 10_CopyNumber_Analysis.md       # CNVkit or GATK CNV
│   └── 11_Clonality_Heterogeneity.md   # PyClone/SciClone
│
├── Libraries_Required.md
├── LICENSE
└── README.md
```
