#!/bin/bash

#SBATCH --job-name=download_sacCer2
#SBATCH --output=download_sacCer2.out
#SBATCH --error=download_sacCer2.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=a.rajalingam.23@abdn.ac.uk

set -euo pipefail

# Output directory
output_dir="/uoa/home/r02ar23/sharedscratch/250113_Hop1_Chip_PRJNA1107023/8_sacCer2_Referance_genome"
mkdir -p "$output_dir"

# Step 1: Download files including md5sum.txt
echo "🔽 Downloading sacCer2 chromosomes and md5sum.txt..."
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/sacCer2/chromosomes/ "${output_dir}/"

cd "$output_dir"

# Step 2: Verify md5 checksums
echo "🔒 Verifying checksums..."
md5sum -c md5sum.txt > md5_check_results.txt || {
    echo "❌ One or more files failed the checksum test. See md5_check_results.txt for details."
    exit 1
}
echo "✅ All files passed MD5 check."

# Step 3: Unzip .fa.gz files
echo "📂 Unzipping FASTA files..."
gunzip -v *.fa.gz

# Step 4: Remove unwanted chromosomes
for chr in 2micron chrM; do
    if [ -f "${chr}.fa" ]; then
        echo "🗑️ Removing ${chr}.fa..."
        rm "${chr}.fa"
    fi
done

# Step 5: Combine in biological order
combined_fasta="sacCer2_combined.fasta"
echo "🛠️ Combining chromosomes into ${combined_fasta} in correct order..."

chrom_order=(chrI chrII chrIII chrIV chrV chrVI chrVII chrVIII chrIX chrX chrXI chrXII chrXIII chrXIV chrXV chrXVI)

# Overwrite if exists
> "${combined_fasta}"

for chr in "${chrom_order[@]}"; do
    if [ -f "${chr}.fa" ]; then
        cat "${chr}.fa" >> "${combined_fasta}"
    else
        echo "❌ Missing chromosome file: ${chr}.fa"
        exit 1
    fi
done

# Step 6: Summary report
echo "📋 Final check: Chromosomes in combined file:"
grep "^>" "${combined_fasta}"

echo "🎯 sacCer2 reference genome preparation completed successfully!"
