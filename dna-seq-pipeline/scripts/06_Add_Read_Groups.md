For Tumor Bulk Sample
```
picard AddOrReplaceReadGroups I= ~/dsa/34/smtl/34sorted.bam O= ~/dsa/34/pcd/34rg.bam RGID=SUB4387364 RGLB=Tumor RGPL=ILLUMINA RGPU=HiSeq2500_1 RGSM=TumorSample1 CREATE_INDEX=true
```

For Normal Blood Sample
```
picard AddOrReplaceReadGroups I= ~/dsa/35/smtl/35sorted.bam O= ~/dsa/35/pcd/35rg.bam RGID=SUB4387365 RGLB=NORMAL RGPL=ILLUMINA RGPU=HiSeq2500_1 RGSM=NormalSample1 CREATE_INDEX=true
```
