```
#For Normal Blood Sample
samtools view -@ 2 -bS ~/dsa/35/bwa/35.sam | \
samtools sort -@ 2 -o ~/dsa/35/smtl/35sorted.bam && \
samtools index -@ 2 -o ~/dsa/35/smtl/35.bai ~/dsa/35/smtl/35sorted.bam

#For Tumor Bulk Sample
samtools view -@ 2 -bS ~/dsa/34/bwa/34.sam | \
samtools sort -@ 2 -o ~/dsa/34/smtl/34sorted.bam && \
samtools index -@ 2 -o ~/dsa/34/smtl/34.bai ~/dsa/34/smtl/34sorted.bam
```
