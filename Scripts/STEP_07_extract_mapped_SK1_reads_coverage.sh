#!/bin/bash -e

#SBATCH --job-name=coverage_mapped_SK1_reads
#SBATCH --output=coverage_mapped_SK1_reads_%j.log
#SBATCH --error=coverage_mapped_SK1_reads_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=a.rajalingam.23@abdn.ac.uk

# Load samtools for coverage calculations
module load samtools

# === Define input and output directories ===
bam_dir="/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/6_mapped_SK1_reads"
coverage_dir="/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/7_mapped_SK1_reads_coverage"
summary_csv="${coverage_dir}/coverage_summary.csv"

# Create output directory if it does not exist
mkdir -p "$coverage_dir"

# === Initialize CSV file to store summary statistics ===
echo "Sample,Total Positions,Covered Positions,Uncovered Positions,Mean Depth,Median Depth" > "$summary_csv"

# === Define list of sample names ===
samples=("IN817" "IN855" "IN819" "IN820" "IP817" "IP855" "IP819" "IP820")

# === Loop through each sample and calculate coverage ===
for sample in "${samples[@]}"; do
    echo "$(date): Calculating coverage for $sample"

    bam_file="${bam_dir}/${sample}_mapped_SK1.bam"
    coverage_file="${coverage_dir}/${sample}_SK1_coverage.txt"
    stats_file="${coverage_dir}/${sample}_SK1_coverage_stats.txt"

    # Skip sample if BAM file does not exist
    if [[ ! -f "$bam_file" ]]; then
        echo "Warning: BAM not found for $sample"
        continue
    fi

    # Generate base-pair level coverage depth
    samtools depth "$bam_file" > "$coverage_file"

    # Compute total, covered, and uncovered positions
    total_positions=$(wc -l < "$coverage_file")
    covered_positions=$(awk '$3 > 0' "$coverage_file" | wc -l)
    uncovered_positions=$((total_positions - covered_positions))

    # Compute mean and median coverage depth
    mean_depth=$(awk '{sum+=$3} END {if (NR>0) print sum/NR; else print 0}' "$coverage_file")
    median_depth=$(awk '{print $3}' "$coverage_file" | sort -n | awk '
        { a[NR] = $1 }
        END {
            if (NR % 2 == 1) print a[(NR + 1) / 2];
            else print (a[NR/2] + a[NR/2 + 1]) / 2;
        }')

    # Write detailed stats to file
    {
        echo "Sample: $sample"
        echo "Total Positions: $total_positions"
        echo "Covered Positions: $covered_positions"
        echo "Uncovered Positions: $uncovered_positions"
        echo "Mean Depth: $mean_depth"
        echo "Median Depth: $median_depth"
    } > "$stats_file"

    # Append stats to summary CSV
    echo "$sample,$total_positions,$covered_positions,$uncovered_positions,$mean_depth,$median_depth" >> "$summary_csv"
    echo "$(date): Done $sample"
done

echo "✅ All SK1-mapped coverage files generated. Summary saved to $summary_csv."
