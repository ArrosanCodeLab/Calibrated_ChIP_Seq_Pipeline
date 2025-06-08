# Step 09: Coverage Calculation of SK1 Reads Mapped to sacCer2

This step involves computing base-level coverage from BAM files of SK1 reads that were previously mapped to the sacCer2 reference genome. Coverage statistics such as total positions, covered/uncovered bases, mean depth, and median depth are generated per sample, and a summary CSV is compiled.

---

## 📁 Required Input and Output Structure

### Input:
Directory of sorted BAM files from step 08:
```
/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/10_sk1_reads_mapped_to_sacCer2/
```

### Output:
Coverage data and statistics will be saved in:
```
/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/11_sk1_reads_mapped_to_sacCer2_coverage/
```

---

## 🧪 Tools Required

- `samtools` (must be loaded via module or available in environment)

---

## 🔁 What This Step Does

For each sample:
1. Load the sorted BAM file aligned to sacCer2.
2. Use `samtools depth -aa` to get base-wise coverage (including positions with 0 reads).
3. Parse coverage data to compute:
   - Total positions
   - Covered positions (coverage > 0)
   - Uncovered positions (coverage = 0)
   - Mean coverage depth
   - Median coverage depth
4. Output sample-specific stats into text files.
5. Append a row for each sample to a summary CSV.

---

## 🖥️ SLURM Script Used

```bash
#!/bin/bash -e

#SBATCH --job-name=calculate_sk1_mapped_to_sacCer2_coverage  # Job name for SLURM
#SBATCH --output=calculate_sk1_mapped_to_sacCer2_coverage_%j.log  # Log for stdout
#SBATCH --error=calculate_sk1_mapped_to_sacCer2_coverage_%j.err   # Log for stderr
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=a.rajalingam.23@abdn.ac.uk

# Load samtools for coverage calculation
module load samtools

# Define input and output directories
mapped_dir="/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/10_sk1_reads_mapped_to_sacCer2"
coverage_dir="/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/11_sk1_reads_mapped_to_sacCer2_coverage"
summary_csv="${coverage_dir}/coverage_summary.csv"

# Create output directory if not exists
mkdir -p "$coverage_dir"

# Write CSV header
echo "Sample,Total Positions,Covered Positions,Uncovered Positions,Mean Depth,Median Depth" > "$summary_csv"

# Define sample names
samples=("IN817" "IN855" "IN819" "IN820" "IP817" "IP855" "IP819" "IP820")

# Loop through each sample
for sample in "${samples[@]}"; do
    echo "$(date): Calculating coverage for sample $sample"

    bam_file="${mapped_dir}/${sample}_mapped_to_sacCer2_sorted.bam"
    coverage_file="${coverage_dir}/${sample}_sk1_mapped_sacCer2_coverage.txt"
    stats_file="${coverage_dir}/${sample}_sk1_mapped_sacCer2_coverage_stats.txt"

    if [[ ! -f "$bam_file" ]]; then
        echo "Error: BAM file for sample $sample not found."
        continue
    fi

    # Generate base-level coverage
    samtools depth -aa "$bam_file" > "$coverage_file"

    # Compute statistics
    total_positions=$(wc -l < "$coverage_file")
    covered_positions=$(awk '$3 > 0' "$coverage_file" | wc -l)
    uncovered_positions=$((total_positions - covered_positions))
    mean_depth=$(awk '{sum += $3} END {printf "%.2f", sum / NR}' "$coverage_file")
    median_depth=$(awk '{print $3}' "$coverage_file" | sort -n | awk '
        { count[NR] = $1 }
        END {
            if (NR % 2 == 1) {
                print count[(NR + 1) / 2]
            } else {
                print (count[NR / 2] + count[(NR / 2 + 1)]) / 2
            }
        }')

    # Save per-sample stats
    {
        echo "Sample: $sample"
        echo "Total Positions: $total_positions"
        echo "Covered Positions: $covered_positions"
        echo "Uncovered Positions: $uncovered_positions"
        echo "Mean Depth: $mean_depth"
        echo "Median Depth: $median_depth"
    } > "$stats_file"

    echo "$sample,$total_positions,$covered_positions,$uncovered_positions,$mean_depth,$median_depth" >> "$summary_csv"
    echo "$(date): Coverage analysis for sample $sample completed."
done

echo "$(date): All coverage calculations completed. Summary saved to $summary_csv."
```

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

## 🧠 Final Note

This final step completes the ChIP-seq processing pipeline by quantifying how well the SK1 reads mapped to the sacCer2 genome. These values are essential for downstream analyses such as peak normalization, differential binding, and reproducibility checks.
