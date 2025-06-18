<h1 align="center">🧬 RNA-Seq Pipeline for Mouse Transcriptomics</h1>

<p align="center">
  <a href="LICENSE">
    <img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License: MIT">
  </a>
</p>

<p align="center">
  <b>Modular, shell-and-R-based RNA-seq analysis workflow for bulk transcriptome analysis in <i>mouse (Mus musculus)</i></b><br>
  <i>Structured for clarity, educational reproducibility, and low-barrier adoption by students or researchers.</i>
</p>

<hr>

<h2>📁 Project Structure</h2>

<pre>
rna-seq-pipeline/
├── config/                   # Configuration files
│   ├── paths.conf
│   ├── samples.csv
│   └── README.md
│
├── data/                     # Raw and processed data (NOT uploaded here)
│   └── README.md             # Instructions to obtain raw data
│
├── scripts/                  # Shell and R scripts for each analysis step
│   ├── 01_Quality-Control.md
│   ├── 02_Trimming.md
│   ├── 03_Genome-Indexing.md
│   ├── 04_Alignment.md
│   ├── 05_SAM-to-BAM.md
│   ├── 06_Sorting.md
│   ├── 07_BAM-Indexing.md
│   ├── 08_featureCount.md
│   └── 09_differential_gene_expression.md
│
├── bash_commands.md          # Optional: useful one-liners
├── Libraries_Required.md     # Required tools and installation
├── LICENSE                   # License file (MIT)
└── README.md                 # Main project documentation
</pre>

<hr>

<h2>🧪 Input Data</h2>

<h3>🔹 Experimental Design</h3>
<ul>
  <li><b>HFD FMT (High-Fat Diet):</b> 2 samples</li>
  <li><b>ND FMT (Normal Diet):</b> 2 samples</li>
  <li>Format: paired-end, mouse fecal microbiota transfer samples.</li>
</ul>

<h3>🔹 Data Sources</h3>
<ul>
  <li>Data downloaded from NCBI SRA using <code>prefetch</code> and <code>fasterq-dump</code>.</li>
  <li>Files are structured like:
    <pre>
data/raw/HFD1_1.fastq.gz
data/raw/HFD1_2.fastq.gz
    </pre>
  </li>
</ul>

<p>
  <b>Sample Configuration:</b><br>
  <a href="https://github.com/VampZie/NGS-Data-Analysis/blob/b3a36de5e152222c6d54e06f289688cdacab4b15/rna-seq-pipeline/config/sample.csv">
    in <code>config/samples.csv</code>
  </a>
</p>

<hr>

<h2>🧬 Reference Files</h2>

<h3>🔸 Reference Genome (FASTA)</h3>
<ul>
  <li><code>GRCm39.primary_assembly.genome.fa</code></li>
  <li>Source: GENCODE Mouse</li>
</ul>

<h3>🔸 Annotation (GFF3)</h3>
<ul>
  <li><code>gencode.vM37.primary_assembly.basic.annotation.gff3</code></li>
  <li>Used in featureCounts for read quantification.</li>
</ul>

<p>
  <a href="https://github.com/VampZie/NGS-Data-Analysis/blob/7c5f0776c35780074f6a0f699f85ba4311f0289b/rna-seq-pipeline/config/paths.conf">
    config/paths.conf
  </a>
</p>

<hr>

<h2>⚙️ Software &amp; Environment</h2>
<p>
  See <a href="https://github.com/VampZie/NGS-Data-Analysis/blob/cb604c83a31ba5f526db246dcb9e25d67d319523/rna-seq-pipeline/Libraries_Required.md"><code>/Libraries_Required.md</code></a> for installation details.
</p>

<hr>

<h2>📊 Analysis Workflow</h2>
<p>
  Each step is documented in <a href="https://github.com/VampZie/NGS-Data-Analysis/tree/a9ab6658ce3ef498704345eb79ba0f903220b381/rna-seq-pipeline/scripts">/script</a> with detailed explanations and commands.
</p>
<ol>
  <li>Quality Control – <code>01_Quality-Control.md</code></li>
  <li>Trimming – <code>02_Trimming.md</code></li>
  <li>Indexing Genome – <code>03_Genome-Indexing.md</code></li>
  <li>Alignment – <code>04_Alignment.md</code></li>
  <li>SAM to BAM – <code>05_SAM-to-BAM.md</code></li>
  <li>Sorting – <code>06_Sorting.md</code></li>
  <li>BAM Indexing – <code>07_BAM-Indexing.md</code></li>
  <li>Quantification – <code>08_featureCount.md</code></li>
  <li>Differential Expression – <code>09_differential_gene_expression.md</code></li>
</ol>

<hr>

<h2>📌 Usage Notes</h2>
<ul>
  <li>🧷 This repository is educational and modular — adapt to your needs.</li>
  <li>🧼 <b>Do not include raw FASTQ or BAM files in your repo.</b> Use <code>.gitignore</code>.</li>
  <li>📜 Scripts are meant to run step-by-step with clarity, not automation (for now).</li>
</ul>

<hr>

<h2>📄 License</h2>
<p>
  This project is licensed under the MIT License – see the <a href="https://github.com/VampZie/NGS-Data-Analysis/blob/f9385362892099133b6d5d70bb84157fb688183b/LICENSE">LICENSE</a> file for details.
</p>

<hr>

<h2>🧑 Author</h2>
<p>
  Maintained by <b>VampZie</b><br>
  <i>Note: This is a personal, not institutional, GitHub identity.</i>
</p>

<hr>

<h2>💬 Feedback</h2>
<p>
  Found a bug, typo, or want to suggest an improvement?<br>
  Please <a href="https://github.com/VampZie/NGS-Data-Analysis/issues">open an issue</a> or <a href="https://github.com/VampZie/NGS-Data-Analysis/discussions">start a discussion</a>.<br>
  📫 You can also reach out via comments on this repository.
</p>

<hr>
