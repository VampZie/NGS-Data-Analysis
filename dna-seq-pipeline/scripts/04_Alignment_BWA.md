```
#For Normal Blood Sample
bwa mem -t 2 ~/dsa/ref/GRCh38.primary_assembly.genome.fa ~/dsa/35/fp/35r1c.fastq ~/dsa/35/fp/35r2c.fastq > ~/dsa/35/bwa/35.sam

#For Tumor Bulk Sample
bwa mem -t 2 ~/dsa/ref/GRCh38.primary_assembly.genome.fa ~/dsa/34/fp/34r1c.fastq ~/dsa/34/fp/34r2c.fastq > ~/dsa/34/bwa/34.sam
```
