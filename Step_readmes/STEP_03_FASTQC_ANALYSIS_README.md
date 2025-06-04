# Step 3: FASTQC Quality Control on Raw and Trimmed FASTQ Files

This step runs FASTQC on both raw and adapter-trimmed paired-end FASTQ files. FASTQC provides quality metrics, helping assess sequence quality before and after trimming.

---

## 🛠️ What This Script Does

- Loads the FASTQC module (ensure it's available on your cluster)
- Runs FASTQC on raw FASTQ files (e.g., `Sample_1.fastq.gz`)
- Runs FASTQC on trimmed FASTQ files (e.g., `Sample_1_trimmed.fastq.gz`)
- Outputs reports in designated directories

---

## 📂 Folder Structure

You should have the following input/output folders:

```
Calibrated_ChIP_Seq/
├── 1_Raw_Sequence_Data/
├── 2_trimmed_output_q20/
└── 3_Fastqc/
    ├── fastqc_raw/
    └── fastqc_trimmed/
```

---

## 🚀 Prerequisites

Ensure FASTQC is available:

```bash
module load fastqc
```

---

## 📜 Script

```bash
#!/bin/bash
#SBATCH --job-name=fastqc_analysis         # Job name
#SBATCH --output=fastqc_analysis.out       # Output file
#SBATCH --error=fastqc_analysis.err        # Error file
#SBATCH --ntasks=1                         # Number of tasks
#SBATCH --cpus-per-task=4                  # Number of CPU cores per task
#SBATCH --mem=32G                          # Memory allocation
#SBATCH --time=1-00:00:00                  # Max runtime (1 day)
#SBATCH --mail-type=BEGIN,END,FAIL         # Email notifications
#SBATCH --mail-user=a.rajalingam.23@abdn.ac.uk  # Notification email

# ======================================================
# Step 3: Run FASTQC on Raw and Trimmed FASTQ Files
# ------------------------------------------------------
# This script performs quality control (QC) on paired-end
# sequencing data using FASTQC. It runs the analysis on
# both raw FASTQ files and adapter-trimmed FASTQ files.
# 
# Input: Raw FASTQ files and trimmed FASTQ files
# Output: FASTQC reports (HTML & ZIP) for each sample
# 
# Sample naming convention:
#   Raw:     <sample_name>_1.fastq.gz, <sample_name>_2.fastq.gz
#   Trimmed: <sample_name>_1_trimmed.fastq.gz, <sample_name>_2_trimmed.fastq.gz
#
# Output directories:
#   - 3_Fastqc/fastqc_raw
#   - 3_Fastqc/fastqc_trimmed
# ======================================================

# Load required FASTQC module
module load fastqc

# Define input and output directories
raw_data_dir="/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/1_Raw_Sequence_Data"
trimmed_data_dir="/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/2_trimmed_output_q20"
fastqc_raw_output_dir="/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/3_Fastqc/fastqc_raw"
fastqc_trimmed_output_dir="/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/3_Fastqc/fastqc_trimmed"

# Create output directories if they don’t exist
mkdir -p "$fastqc_raw_output_dir"
mkdir -p "$fastqc_trimmed_output_dir"

# Define sample names (update if your samples differ)
sample_names=("IN817" "IN855" "IN819" "IN820" "IP817" "IP855" "IP819" "IP820")

# Run FASTQC on raw FASTQ files
echo "Running FASTQC on raw FASTQ files..."
for i in "${sample_names[@]}"
do
    echo "▶ Processing raw data: ${i}"
    fastqc "${raw_data_dir}/${i}_1.fastq.gz" \
           "${raw_data_dir}/${i}_2.fastq.gz" \
           -o "$fastqc_raw_output_dir"
done

# Run FASTQC on trimmed FASTQ files
echo "Running FASTQC on trimmed FASTQ files..."
for i in "${sample_names[@]}"
do
    echo "▶ Processing trimmed data: ${i}"
    fastqc "${trimmed_data_dir}/${i}_1_trimmed.fastq.gz" \
           "${trimmed_data_dir}/${i}_2_trimmed.fastq.gz" \
           -o "$fastqc_trimmed_output_dir"
done

echo "✅ FASTQC analysis completed for both raw and trimmed data."
```

---

## 👨‍🔬 Author
**Arrosan Rajalingam**  
PhD Student – Murakami Lab  
University of Aberdeen, UK  
📧 a.rajalingam.23@abdn.ac.uk

## 🤝 Contributor
**Yusuke Tsuruta**  
Tokyo Metropolitan University, Japan

## 📬 Corresponding Author
**Hajime Murakami**  
University of Aberdeen, UK  
📧 hajime.murakami1@abdn.ac.uk
