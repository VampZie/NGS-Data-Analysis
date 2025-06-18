<p align="center">
  <img src="https://img.icons8.com/fluency/96/dna-helix.png" width="90" alt="RNA-Seq">
</p>

<h1 align="center" style="font-size:2.6em; font-weight:bold;">🧬 RNA-Seq Pipeline for Mouse Transcriptomics</h1>

<p align="center">
  <a href="LICENSE">
    <img src="https://img.shields.io/badge/License-MIT-yellow.svg" alt="License: MIT">
  </a>
</p>
<p align="center" style="font-size:1.2em;">
  <b>Modular, shell-and-R-based RNA-seq analysis workflow for bulk transcriptome analysis in <i>Mus musculus</i></b><br>
  <i>Structured for clarity, reproducibility, and easy adoption by students & researchers.</i>
</p>

<hr style="height:2px; border:none; background:linear-gradient(to right,#6a82fb,#fc5c7d);">

<h2>📦 Project Structure</h2>

<pre style="background-color:#f8f8f8; border-radius:7px; padding:16px; border: 1px solid #e3e3e3;">
rna-seq-pipeline/
├── <b>config/</b>                   # Configuration files
│   ├── paths.conf
│   ├── samples.csv
│   └── README.md
│
├── <b>data/</b>                     # Raw and processed data (NOT uploaded here)
│   └── README.md             # Instructions to obtain raw data
│
├── <b>scripts/</b>                  # Shell and R scripts for each analysis step
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

---

<h2>🧪 Input Data</h2>

<table>
  <tr>
    <td width="50%">
      <h3>🔹 <b>Experimental Design</b></h3>
      <ul>
        <li><b>HFD FMT (High-Fat Diet):</b> 2 samples</li>
        <li><b>ND FMT (Normal Diet):</b> 2 samples</li>
        <li>Format: <b>paired-end</b>, mouse fecal microbiota transfer samples</li>
      </ul>
    </td>
    <td width="50%">
      <h3>🔹 <b>Data Sources</b></h3>
      <ul>
        <li>Download from NCBI SRA via <code>prefetch</code> & <code>fasterq-dump</code></li>
        <li>File pattern:
          <br><code>data/raw/HFD1_1.fastq.gz<br>data/raw/HFD1_2.fastq.gz</code>
        </li>
      </ul>
    </td>
  </tr>
</table>

<blockquote>
  <b>Sample configuration:</b> See
  <a href="https://github.com/VampZie/NGS-Data-Analysis/blob/b3a36de5e152222c6d54e06f289688cdacab4b15/rna-seq-pipeline/config/sample.csv">
    <code>config/samples.csv</code>
  </a>
</blockquote>

---

<h2>🧬 Reference Files</h2>

<ul>
  <li>
    <b>Reference Genome (FASTA):</b> <code>GRCm39.primary_assembly.genome.fa</code>
    <br><small>Source: GENCODE Mouse</small>
  </li>
  <li>
    <b>Annotation (GFF3):</b> <code>gencode.vM37.primary_assembly.basic.annotation.gff3</code>
    <br><small>Used in featureCounts for read quantification.</small>
  </li>
</ul>
<p>
  <a href="https://github.com/VampZie/NGS-Data-Analysis/blob/7c5f0776c35780074f6a0f699f85ba4311f0289b/rna-seq-pipeline/config/paths.conf">
    <code>config/paths.conf</code>
  </a>
</p>

---

<h2>⚙️ Software & Environment</h2>
<p>
  See <a href="https://github.com/VampZie/NGS-Data-Analysis/blob/cb604c83a31ba5f526db246dcb9e25d67d319523/rna-seq-pipeline/Libraries_Required.md"><code>Libraries_Required.md</code></a> for installation details.
</p>

---

<h2>📊 Analysis Workflow</h2>
<p>
  <img src="https://img.icons8.com/color/48/workflow.png" width="28" align="top"/> Each step is documented in <a href="https://github.com/VampZie/NGS-Data-Analysis/tree/a9ab6658ce3ef498704345eb79ba0f903220b381/rna-seq-pipeline/scripts"><code>scripts/</code></a> with explanations & commands.
</p>
<ol>
  <li><b>Quality Control</b> – <code>01_Quality-Control.md</code></li>
  <li><b>Trimming</b> – <code>02_Trimming.md</code></li>
  <li><b>Indexing Genome</b> – <code>03_Genome-Indexing.md</code></li>
  <li><b>Alignment</b> – <code>04_Alignment.md</code></li>
  <li><b>SAM to BAM</b> – <code>05_SAM-to-BAM.md</code></li>
  <li><b>Sorting</b> – <code>06_Sorting.md</code></li>
  <li><b>BAM Indexing</b> – <code>07_BAM-Indexing.md</code></li>
  <li><b>Quantification</b> – <code>08_featureCount.md</code></li>
  <li><b>Differential Expression</b> – <code>09_differential_gene_expression.md</code></li>
</ol>

---

<h2>📌 Usage Notes</h2>
<ul>
  <li>🧷 <b>Educational & modular</b> — adapt to your needs.</li>
  <li>🧼 <b>Do NOT commit raw FASTQ or BAM files</b> — use <code>.gitignore</code>.</li>
  <li>📜 Scripts are designed to run step-by-step for clarity (not automation, for now).</li>
</ul>

---

<h2>📄 License</h2>
<p>
  <img src="https://img.icons8.com/color/32/license.png" width="20"/> MIT License – see
  <a href="https://github.com/VampZie/NGS-Data-Analysis/blob/f9385362892099133b6d5d70bb84157fb688183b/LICENSE"><b>LICENSE</b></a>
</p>

---

<h2>🧑 Author</h2>
<table>
  <tr>
    <td width="35">
      <img src="https://img.icons8.com/color/48/cat-profile.png" width="32"/>
    </td>
    <td>
      <b>Maintained by:</b> VampZie<br>
      <i>Note: This is a personal, not institutional, GitHub identity.</i>
    </td>
  </tr>
</table>

---

<h2>💬 Feedback</h2>
<p>
  <img src="https://img.icons8.com/fluency/32/feedback.png" width="20"/> Found a bug, typo, or want to suggest an improvement?<br>
  👉 <a href="https://github.com/VampZie/NGS-Data-Analysis/issues">Open an issue</a>
  or <a href="https://github.com/VampZie/NGS-Data-Analysis/discussions">start a discussion</a>.<br>
  📫 Or reach out via repository comments.
</p>

<hr style="height:2px; border:none; background:linear-gradient(to right,#fc5c7d,#6a82fb);">
