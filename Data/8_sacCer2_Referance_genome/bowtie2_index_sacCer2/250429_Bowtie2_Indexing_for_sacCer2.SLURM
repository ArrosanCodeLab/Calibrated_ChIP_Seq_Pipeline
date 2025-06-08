#!/bin/bash

#SBATCH --job-name=bowtie2_index_sacCer2
#SBATCH --output=bowtie2_index_sacCer2.out
#SBATCH --error=bowtie2_index_sacCer2.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=a.rajalingam.23@abdn.ac.uk

set -e
set -u
set -o pipefail

# Paths
fasta_file="/uoa/home/r02ar23/sharedscratch/250113_Hop1_Chip_PRJNA1107023/8_sacCer2_Referance_genome/sacCer2_combined.fasta"
index_dir="/uoa/home/r02ar23/sharedscratch/250113_Hop1_Chip_PRJNA1107023/8_sacCer2_Referance_genome/bowtie2_index_sacCer2"
index_prefix="${index_dir}/sacCer2_combined"

# Create index directory if it doesn't exist
mkdir -p "${index_dir}"

# Build Bowtie2 index
echo "🔨 Building Bowtie2 index for sacCer2_combined.fasta..."
bowtie2-build -f "${fasta_file}" "${index_prefix}"

# Check if index built successfully
if [ $? -eq 0 ]; then
    echo "✅ Bowtie2 index created successfully at: ${index_prefix}.*"
else
    echo "❌ Bowtie2 indexing failed!"
    exit 1
fi

echo "🎯 Script finished!"
