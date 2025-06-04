#!/bin/zsh

# ===============================================================
# 🧬 Script: build_bowtie2_index.zsh
# 🔧 Purpose: To build a Bowtie2 index from the combined SK1 + Smik genome
# 📦 Required: Bowtie2 installed and accessible via command line
# 🔄 Input: Combined genome FASTA file (SK1 + Smik)
# 🧾 Output: Bowtie2 index files (.bt2) in a dedicated directory
# ===============================================================

# ======================= IMPORTANT ==============================
# 🔁 Before using this script:
# 1. Ensure that you have downloaded the folder `3_reference_genomes`
#    from the GitHub repository and placed it in:
#    /uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq
# 2. That folder must contain:
#    - genome_seq_SK1_and_Smik.fa
#    - SK1/ and Smik/ folders with original references (optional)
# ==============================================================

# ===> Path to the combined genome FASTA file <===
genome_fasta="/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/4_reference_genomes/genome_seq_SK1_and_Smik.fa"

# ===> Output directory where Bowtie2 index will be stored <===
output_index="/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/4_reference_genomes/combined_genome_index"

# ===> Create the output directory if it doesn’t already exist <===
mkdir -p "${output_index}"

# ===> Run Bowtie2 to build index <===
echo "🔨 Building Bowtie2 index for ${genome_fasta}..."
bowtie2-build -f "${genome_fasta}" "${output_index}/genome_seq_SK1_and_Smik"

# ===> Check for successful completion and notify user <===
if [ $? -eq 0 ]; then
    echo "✅ Bowtie2 index building completed successfully."
else
    echo "❌ Error: Bowtie2 index building failed." >&2
fi
