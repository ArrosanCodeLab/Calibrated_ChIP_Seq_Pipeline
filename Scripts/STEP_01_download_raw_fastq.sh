#!/bin/bash
#SBATCH --job-name=download_fastq
#SBATCH --output=download_fastq_%j.out
#SBATCH --error=download_fastq_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=a.rajalingam.23@abdn.ac.uk

# ======================= IMPORTANT =======================
# Update only the two lines marked “==> UPDATE THIS <==”.
# =========================================================

# Define the directory for raw data downloads
# ===> we want to create 1_Raw_Sequence_Data under the main folder <===
raw_data_dir="/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/1_Raw_Sequence_Data"

# Define the path to the metadata file (MD5 + URLs)
# ===> metadata is located one level up, in the main Calibrated_ChIP_Seq folder <===
md5sum_file="/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/EN****_22samples_md5sum_DownloadLink.txt"

# ---------------------------------------------------------
# Do NOT modify anything below unless you know what you’re doing.
# ---------------------------------------------------------

# Ensure the raw data directory exists (will create 
# /uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/1_Raw_Sequence_Data)
mkdir -p "$raw_data_dir"

# Change to the raw data directory
cd "$raw_data_dir" || { echo "Cannot change to raw data directory!"; exit 1; }

# Download each file listed in the metadata file
echo "Starting downloads..."
awk 'NR>1 {print $4}' "$md5sum_file" | while read -r url; do
    filename=$(basename "$url")
    if [ -f "$filename" ]; then
        echo "File $filename already exists. Skipping download."
    else
        echo "Downloading $url..."
        wget "$url" || { echo "Failed to download $url"; exit 1; }
    fi
done

# Validate downloads using MD5 checksums
echo "Validating downloads with MD5 checksums..."
awk 'NR>1 {print $3, $1}' "$md5sum_file" | while read -r checksum filename; do
    if md5sum -c <(echo "$checksum $filename"); then
        echo "$filename: VALID"
    else
        echo "$filename: INVALID"
        exit 1
    fi
done

echo "All files downloaded and validated successfully!"
