# 🧬 Calibrated ChIP-Seq Pipeline – Beginner-Friendly Guide

This repository contains a **step-by-step pipeline** for analyzing *calibrated ChIP-seq data* in yeast. It is designed to be **modular**, **reproducible**, and **accessible** even if you’re new to ChIP-seq or bioinformatics.

---

## 🧠 What Is ChIP-Seq?

**ChIP-seq (Chromatin Immunoprecipitation Sequencing)** is a method to find where specific proteins bind on DNA. It’s widely used to study:
- Histone modifications
- Transcription factors
- Structural proteins like cohesin (e.g., Rec8)

---

## 📏 What Is Calibrated ChIP-Seq?

Standard ChIP-seq can suffer from technical variability (e.g., differences in sequencing depth). **Calibrated ChIP-seq** improves accuracy by adding **spike-in cells** from a different species.

In this pipeline:
- *Saccharomyces cerevisiae* (SK1) is the main organism.
- *Saccharomyces mikatae* (Smik) is the spike-in control.

---

## 🧪 What Does This Pipeline Do?

This pipeline takes raw sequencing reads and produces **normalized binding profiles** across the genome.

### 🔬 Main Steps:
1. Download raw data and check file integrity.
2. Trim sequencing adapters and low-quality reads.
3. Run quality control (FastQC).
4. Map reads to a combined genome (SK1 + Smik).
5. Separate reads mapped to SK1 only.
6. Calculate genome-wide coverage.
7. Remap SK1 reads to a standard yeast genome (sacCer2).
8. Generate normalized coverage statistics.

---

## 📊 Pipeline Flowchart

```
            +------------------+
            |  Raw FASTQ Files |
            +------------------+
                     ↓
              Adapter Trimming
                     ↓
               Quality Control
                     ↓
        Mapping to SK1 + Smik Genome
                     ↓
         ┌────────────┴────────────┐
         ↓                         ↓
  Reads Mapped to SK1       Reads Mapped to Smik
         ↓
   Remap to sacCer2 Genome
         ↓
   Final Coverage and Stats
```

---

## 📁 Repository Structure

```
Calibrated_ChIP_Seq_Pipeline/
├── Data/
│   ├── 1_Raw_Sequence_Data/                      # Raw input FASTQs
│   ├── 2_trimmed_output_q20/                     # Adapter-trimmed reads
│   ├── 3_Fastqc/
│   │   ├── fastqc_raw/                           # FastQC reports (raw)
│   │   └── fastqc_trimmed/                       # FastQC reports (trimmed)
│   ├── 4_reference_genomes/                      # Combined genome + index
│   ├── 5_Combain_genome_mapping_output/          # Bowtie2 mapped BAMs + stats
│   ├── 6_mapped_SK1_reads/                       # SK1-only filtered BAMs/FASTQs
│   ├── 7_mapped_SK1_reads_coverage/              # Base-level coverage from SK1 BAMs
│   ├── 8_sacCer2_Referance_genome/               # sacCer2 FASTA and Bowtie2 index
│   ├── 9_sk1_reads_mapped_to_sacCer2/            # Final SK1 reads aligned to sacCer2
│   └── 10_sk1_reads_mapped_to_sacCer2_coverage/  # Final normalized coverage data
├── Scripts/                                      # All SLURM-compatible shell scripts
├── Step_readmes/                                 # Individual README for each step
└── README.md                                     # Top-level pipeline summary
```

---

## 🚀 How to Use the Pipeline

1. **Clone the repository**
```bash
git clone https://github.com/ArrosanCodeLab/Calibrated_ChIP_Seq_Pipeline.git
cd Calibrated_ChIP_Seq_Pipeline
```

2. **Place your raw FASTQ files** in `Data/1_Raw_Sequence_Data/`.

3. **Make sure all reference genomes** are in the right folders under `Data/`.

4. **Run each step** using SLURM:
```bash
sbatch Scripts/STEP_01_download_raw_fastq.sh
sbatch Scripts/STEP_02_adapter_trimming.sh
...
```

5. **Check `Step_readmes/`** for full explanation and output examples for each step.

---

## 🔁 Full Workflow Steps

Each step has an associated SLURM script and README.

| Step | Description | Script |
|------|-------------|--------|
| 01   | Download raw FASTQ files & validate MD5 | `STEP_01_download_raw_fastq.sh` |
| 02   | Adapter trimming using `fastp` | `STEP_02_adapter_trimming.sh` |
| 03   | Quality control using `FastQC` | `STEP_03_FastQC_analysis.sh` |
| 04   | Build Bowtie2 index for SK1+Smik | `STEP_04_build_bowtie2_index.zsh` |
| 05   | Map reads to combined genome | `STEP_05_bowtie2_mapping_combain_referance_genome.sh` |
| 05A  | Generate detailed mapping summary | `STEP_05A_bowtie2_detailed_mapping_summary.sh` |
| 05B  | Tabular mapping summary | `STEP_05B_bowtie2_detailed_mapping_summary_table.sh` |
| 06   | Extract SK1-only reads | `STEP_06_extract_mapped_SK1_read_as_fastq.sh` |
| 07   | Compute SK1 read coverage | `STEP_07_extract_mapped_SK1_reads_coverage.sh` |
| 08   | Remap SK1 reads to sacCer2 | `STEP_08_sk1_reads_mapped_to_sacCer2.sh` |
| 09   | Final sacCer2 coverage | `STEP_09_Sk1_reads_mapped_sacCer2_Coverage.sh` |

---

## 🛠 Requirements
- Linux or Mac OS
- SLURM-based HPC system
- Modules:
  - `fastp`
  - `fastqc`
  - `bowtie2`
  - `samtools`
  - `bedtools2`

---

## 💡 Who Is This For?

This pipeline is ideal for:
- Students new to ChIP-seq
- Wet-lab biologists analyzing ChIP data for the first time
- Yeast researchers using SK1 with Smik spike-ins

---

## 👨‍🔬 Author

**Arrosan Rajalingam**  
PhD Student – Murakami Lab  
University of Aberdeen, UK  
📧 a.rajalingam.23@abdn.ac.uk

---

## 🤝 Contributor

**Yusuke Tsuruta**  
Tokyo Metropolitan University, Japan

---

## 📬 Corresponding Author

**Dr. Hajime Murakami**  
University of Aberdeen, UK  
📧 hajime.murakami1@abdn.ac.uk

---

## ⚠️ Note

We’re actively improving this pipeline. If you face any issues or need help, please open an issue on GitHub or email the maintainer. An automated version of this pipeline is under development!
