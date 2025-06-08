# Step 7: Generate Coverage from SK1-Mapped Reads

This step calculates genome-wide coverage from BAM files that contain reads mapped specifically to SK1 chromosomes. It uses `samtools depth` to compute depth per position and generates both per-sample coverage files and a summary table.

---

## 📌 What This Step Does

* Loads the `samtools` module
* Iterates over all SK1-specific BAM files
* Computes depth using `samtools depth`
* Calculates statistics: total positions, covered positions, uncovered positions, mean and median depth
* Outputs:

  * A detailed `.txt` file with position-wise coverage
  * A `.txt` file with summary statistics
  * A combined CSV summary across all samples

---

## 📁 Input & Output Structure

* **Input Directory:**

  ```
  /uoa/home/r02ar23/sharedscratch/250502_Rec8_CTD_Red1_helix_Rec8-6myc_ChIP_Sq/6_250503_mapped_SK1_reads
  ```

  * Contains `${sample}_mapped_SK1.bam` files

* **Output Directory:**

  ```
  /uoa/home/r02ar23/sharedscratch/250429_Rec8_Chimeric_Construct_V5-Red1_ChIP-Sq/7_250503_mapped_SK1_reads_coverage
  ```

  * Will contain:

    * Coverage depth files: `${sample}_SK1_coverage.txt`
    * Stats summary: `${sample}_SK1_coverage_stats.txt`
    * Combined summary: `coverage_summary.csv`

---

## 🧪 Tools Required

* `samtools` (should be loaded via module)

---

## 📜 Full Shell Script with Descriptions

```bash
#!/bin/bash -e

#SBATCH --job-name=coverage_mapped_SK1_reads  # Set SLURM job name
#SBATCH --output=coverage_mapped_SK1_reads_%j.log  # SLURM stdout
#SBATCH --error=coverage_mapped_SK1_reads_%j.err   # SLURM stderr
#SBATCH --ntasks=1                      # Number of tasks (1 for serial)
#SBATCH --cpus-per-task=4               # Number of CPU cores
#SBATCH --mem=16G                       # Memory allocation
#SBATCH --time=1-00:00:00               # Max runtime (1 day)
#SBATCH --mail-type=BEGIN,END,FAIL      # Notifications
#SBATCH --mail-user=a.rajalingam.23@abdn.ac.uk  # Email for notifications

# Load samtools module
module load samtools

# Input: BAM files mapped only to SK1
bam_dir="/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/6_mapped_SK1_reads"

# Output directory for coverage data
coverage_dir="/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/7_mapped_SK1_reads_coverage"

# Summary file to collect per-sample stats
summary_csv="${coverage_dir}/coverage_summary.csv"

# Create output directory if it doesn't exist
mkdir -p "$coverage_dir"

# Initialize CSV header
echo "Sample,Total Positions,Covered Positions,Uncovered Positions,Mean Depth,Median Depth" > "$summary_csv"

# Sample list (adjust as needed)
samples=("IN817" "IN855" "IN819" "IN820" "IP817" "IP855" "IP819" "IP820")

# Loop over each sample
for sample in "${samples[@]}"; do
    echo "$(date): Calculating coverage for $sample"

    bam_file="${bam_dir}/${sample}_mapped_SK1.bam"
    coverage_file="${coverage_dir}/${sample}_SK1_coverage.txt"
    stats_file="${coverage_dir}/${sample}_SK1_coverage_stats.txt"

    # Check if BAM file exists
    if [[ ! -f "$bam_file" ]]; then
        echo "Warning: BAM not found for $sample"
        continue
    fi

    # Generate position-wise coverage using samtools depth
    samtools depth "$bam_file" > "$coverage_file"

    # Calculate total covered positions, mean, median etc.
    total_positions=$(wc -l < "$coverage_file")
    covered_positions=$(awk '$3 > 0' "$coverage_file" | wc -l)
    uncovered_positions=$((total_positions - covered_positions))
    mean_depth=$(awk '{sum+=$3} END {if (NR>0) print sum/NR; else print 0}' "$coverage_file")
    median_depth=$(awk '{print $3}' "$coverage_file" | sort -n | awk '
        { a[NR] = $1 }
        END {
            if (NR % 2 == 1) print a[(NR + 1) / 2];
            else print (a[NR/2] + a[NR/2 + 1]) / 2;
        }')

    # Save individual stats
    {
        echo "Sample: $sample"
        echo "Total Positions: $total_positions"
        echo "Covered Positions: $covered_positions"
        echo "Uncovered Positions: $uncovered_positions"
        echo "Mean Depth: $mean_depth"
        echo "Median Depth: $median_depth"
    } > "$stats_file"

    # Append to summary CSV
    echo "$sample,$total_positions,$covered_positions,$uncovered_positions,$mean_depth,$median_depth" >> "$summary_csv"

    echo "$(date): Done $sample"
done

echo "✅ All SK1-mapped coverage files generated. Summary saved to $summary_csv."
```

---

## 👨‍🔬 Author

**Arrosan Rajalingam**
PhD Student – Murakami Lab
University of Aberdeen, UK
📧 [a.rajalingam.23@abdn.ac.uk](mailto:a.rajalingam.23@abdn.ac.uk)

## 🤝 Contributor

**Yusuke Tsuruta**
Tokyo Metropolitan University, Japan

## 📬 Corresponding Author

**Dr. Hajime Murakami**
University of Aberdeen, UK
📧 [hajime.murakami1@abdn.ac.uk](mailto:hajime.murakami1@abdn.ac.uk)
