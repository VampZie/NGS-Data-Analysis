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
│   ├── 01_Quality_Control.md
│   ├── 02_Trimming.md
│   ├── 03_Reference_Genome_Indexing.md
│   ├── 04_Alignment_BWA.md
│   ├── 05_Sort_Index_BAM.md
│   ├── 06_Add_Read_Groups.md        
│   ├── 07_Mark_Duplicates.md
│   ├── 08_BaseRecalibration.md
│   ├── 09_Somatic_Variant_Calling.md   # GATK Mutect2
│   ├── 10_Filtering_Variants.md
│   ├── 11_Annotation.md                # VEP or ANNOVAR
│   ├── 12_CopyNumber_Analysis.md       # CNVkit or GATK CNV
│   ├── 13_Clonality_Heterogeneity.md   # PyClone/SciClone
│
├── Libraries_Required.md
├── LICENSE
└── README.md
```
