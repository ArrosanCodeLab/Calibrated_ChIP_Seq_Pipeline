#!/bin/bash -e

#SBATCH --job-name=calculate_sk1_mapped_to_sacCer2_coverage
#SBATCH --output=calculate_sk1_mapped_to_sacCer2_coverage_%j.log
#SBATCH --error=calculate_sk1_mapped_to_sacCer2_coverage_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8 
#SBATCH --mem=128G 
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=a.rajalingam.23@abdn.ac.uk

# Load required modules
module load samtools

# Define directories
mapped_dir="/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/9_sk1_reads_mapped_to_sacCer2"
coverage_dir="/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/10_sk1_reads_mapped_to_sacCer2_coverage"
summary_csv="${coverage_dir}/coverage_summary.csv"

# Create coverage output directory
mkdir -p "$coverage_dir"

# Initialize summary CSV
echo "Sample,Total Positions,Covered Positions,Uncovered Positions,Mean Depth,Median Depth" > "$summary_csv"

# Define actual sample names
samples=("IN817" "IN855" "IN819" "IN820" "IP817" "IP855" "IP819" "IP820")

# Loop through samples
for sample in "${samples[@]}"; do
    echo "$(date): Calculating coverage for sample $sample"

    # Define input and output files
    bam_file="${mapped_dir}/${sample}_mapped_to_sacCer2_sorted.bam"
    coverage_file="${coverage_dir}/${sample}_sk1_mapped_sacCer2_coverage.txt"
    stats_file="${coverage_dir}/${sample}_sk1_mapped_sacCer2_coverage_stats.txt"

    # Check if BAM file exists
    if [[ ! -f "$bam_file" ]]; then
        echo "Error: BAM file for sample $sample not found: $bam_file"
        continue
    fi

    # Calculate coverage using samtools depth
    samtools depth -aa "$bam_file" > "$coverage_file"

    # Calculate coverage statistics
    echo "$(date): Calculating coverage statistics for sample $sample"
    total_positions=$(wc -l < "$coverage_file")
    covered_positions=$(awk '$3 > 0' "$coverage_file" | wc -l)
    uncovered_positions=$((total_positions - covered_positions))
    mean_depth=$(awk '{sum += $3} END {printf "%.2f", sum / NR}' "$coverage_file")

    # Calculate median depth
    median_depth=$(awk '{print $3}' "$coverage_file" | sort -n | awk '
    {
        count[NR] = $1;
    }
    END {
        if (NR % 2 == 1) {
            median = count[(NR + 1) / 2];
        } else {
            median = (count[NR / 2] + count[(NR / 2) + 1]) / 2;
        }
        printf "%.2f", median;
    }')

    # Write sample-specific stats to file
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

    echo "$(date): Coverage analysis for sample $sample completed."
done

# Completion message
echo "$(date): Coverage analysis for all samples completed."
echo "Summary saved to $summary_csv."
