for BaseRecallibration you need two files one is reference.fa.fai and dict
```
#For creating the dict
picard CreateSequenceDictionary   R=~/dsa/ref/GRCh38.primary_assembly.genome.fa   O=~/dsa/ref/GRCh38.primary_assembly.genome.dict

#For creating the fa.fai
samtools faidx GRCh38.primary_assembly.genome.fa

```
