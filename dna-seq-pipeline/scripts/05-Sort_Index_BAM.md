```
samtools view -@ 2 -bS ~/dsa/35/bwa/35.sam | \
samtools sort -@ 2 -o ~/dsa/35/smtl/35sorted.bam && \
samtools index -@ 2 -o ~/dsa/35/smtl/35.bai ~/dsa/35/smtl/35sorted.bam
```
