#!/bin/bash
#SBATCH --job-name=extract_mapped_SK1     # Job name
#SBATCH --output=extract_mapped_SK1.out   # Standard output log
#SBATCH --error=extract_mapped_SK1.err    # Standard error log
#SBATCH --ntasks=1                        # Number of tasks
#SBATCH --cpus-per-task=8                 # CPUs per task (increased from 4 to 8)
#SBATCH --mem=128G                        # Memory allocation (increased from 16G to 128G)
#SBATCH --time=1-00:00:00                 # Maximum runtime (1 day)
#SBATCH --mail-type=BEGIN,END,FAIL        # Email notifications
#SBATCH --mail-user=a.rajalingam.23@abdn.ac.uk  # User email for notifications

# Define paths
input_dir="/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/5_Combain_genome_mapping_output"
output_dir="/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/6_mapped_SK1_reads"

# Create output directory if not exists
mkdir -p "${output_dir}"

# Array of sample names matching the combined genome mapping script
sample_names=("IN817" "IN855" "IN819" "IN820" "IP817" "IP855" "IP819" "IP820")

# Extract reads mapped to SK1 chromosomes
for sample in "${sample_names[@]}"; do
    echo "Processing sample: ${sample}"

    input_bam="${input_dir}/${sample}_mapped.bam"
    output_bam="${output_dir}/${sample}_mapped_SK1.bam"
    output_fastq="${output_dir}/${sample}_mapped_SK1.fq.gz"

    # Filter reads mapped to SK1 chromosomes only
    echo "Filtering reads mapped to SK1 chromosomes..."
    samtools view -@ 7 -b "${input_bam}" SK1_1 SK1_2 SK1_3 SK1_4 SK1_5 SK1_6 SK1_7 SK1_8 SK1_9 SK1_10 SK1_11 SK1_12 SK1_13 SK1_14 SK1_15 SK1_16 > "${output_bam}"

    # Index the filtered BAM file
    echo "Indexing SK1-specific BAM file..."
    samtools index "${output_bam}"

    # Convert SK1-specific BAM to FASTQ and compress it
    echo "Generating FASTQ file..."
    samtools fastq -@ 7 "${output_bam}" | gzip > "${output_fastq}"

    echo "Completed processing for sample: ${sample}"
done

echo "All samples processed successfully!"
