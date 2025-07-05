For Tumor Bulk Sample
```
picard MarkDuplicates I=/home/vzscyborg/dsa/34/pcd/34rg.bam O=/home/vzscyborg/dsa/34/pcd/34dup.bam M=/home/vzscyborg/dsa/34/pcd/34dup_mtrx.txt CREATE_INDEX=true
```

For Normal Blood Sample
```
picard MarkDuplicates I=/home/vzscyborg/dsa/35/pcd/35rg.bam O=/home/vzscyborg/dsa/35/pcd/35dup.bam M=/home/vzscyborg/dsa/35/pcd/35dup_mtrx.txt CREATE_INDEX=true
```


As this step already generated ```indexed bam``` file so i don't need to do indexing again if you wont find ```'xxxx.bai' and ```'xxxx.bam'``` file, do samtool indexing command.
