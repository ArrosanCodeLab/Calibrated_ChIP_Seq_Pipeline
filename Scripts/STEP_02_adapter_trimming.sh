#!/bin/bash
#SBATCH --job-name=adapter_trimming             # Job name
#SBATCH --output=adapter_trimming_%j.out        # Output file
#SBATCH --error=adapter_trimming_%j.err         # Error file
#SBATCH --ntasks=1                              # Number of tasks
#SBATCH --cpus-per-task=4                       # Number of CPU cores per task
#SBATCH --mem=16G                               # Memory allocation
#SBATCH --time=04:00:00                         # Max runtime
#SBATCH --mail-type=BEGIN,END,FAIL              # Email notifications for job status
#SBATCH --mail-user=a.rajalingam.23@abdn.ac.uk  # Your email for job notifications

# ======================= IMPORTANT =======================
# Update only the two sections marked “==> UPDATE THIS <==”.
# =========================================================

# Define the output directory for trimmed files
# ===> UPDATE THIS to your desired trimmed output directory <===
output_dir="/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/2_trimmed_output_q20"
mkdir -p "$output_dir"  # Ensure the directory exists

# Define your sample names for paired-end reads
# ===> UPDATE THIS to reflect your actual sample names <===
sample_names=("IN817" "IN855" "IN819" "IN820" "IP817" "IP855" "IP819" "IP820")

# ---------------------------------------------------------
# Do NOT modify anything below unless you know what you're doing.
# ---------------------------------------------------------

# Load the required module for fastp
module load fastp

# Loop over all samples and trim adapters using fastp
# Run fastp for each paired-end sample
    # ===> Make sure your FASTQ filenames follow this pattern: ${i}_1.fastq.gz and ${i}_2.fastq.gz
    # ===> If your filenames differ, UPDATE the suffix accordingly (e.g., _R1.fastq.gz or .fq.gz)
for i in "${sample_names[@]}"
do
    echo "Processing sample: ${i}"

    fastp -i "/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/1_Raw_Sequence_Data/${i}_1.fastq.gz" \
          -I "/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/1_Raw_Sequence_Data/${i}_2.fastq.gz" \
          -3 \
          -o "$output_dir/${i}_1_trimmed.fastq.gz" \
          -O "$output_dir/${i}_2_trimmed.fastq.gz" \
          -q 20 -l 20

    echo "Completed processing for sample: ${i}"
done

echo "Adapter trimming completed for all samples."
