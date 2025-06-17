## 🔧 Setup Instructions

Before running the pipeline:

1. Download the following files from [GENCODE](https://www.gencodegenes.org/mouse/):
   - Genome FASTA: `GRCm39.primary_assembly.genome.fa`
   - Annotation GFF3: `gencode.vM37.primary_assembly.basic.annotation.gff3`

2. Place them in the following directory:
   - data/reference/

3. Confirm your `config/paths.conf` contains:
```bash
GENOME=data/reference/GRCm39.primary_assembly.genome.fa
ANNOTATION=data/reference/gencode.vM37.primary_assembly.basic.annotation.gff3
```
