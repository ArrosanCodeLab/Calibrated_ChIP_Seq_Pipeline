# Step 08: Map SK1 Reads to sacCer2 Reference Genome

This step involves taking the SK1-specific reads (already extracted in a previous step) and remapping them to the `sacCer2` reference genome using Bowtie2. This is important to analyze only those reads that belong to the SK1 strain in the context of the standard yeast genome annotation.

---

## 📁 Required Directory Structure

Make sure the following folder is present **in your main working directory**:

```
8_sacCer2_Referance_genome/
```

This folder **must contain**:
- All sacCer2 individual chromosome FASTA files (e.g., `chrI.fa`, `chrII.fa`, ..., `chrXVI.fa`)
- A combined reference genome file `sacCer2_combined.fasta`
- Bowtie2 index files in a subdirectory: `bowtie2_index_sacCer2/`
- Validation and log files (e.g., `md5sum.txt`, `README.txt`, etc.)

---

## 📌 What This Step Does

For each SK1 sample:
1. Sort BAM files by read name (for paired FASTQ conversion).
2. Convert sorted BAM to paired FASTQ (`R1.fastq` and `R2.fastq`).
3. Validate the number of reads in R1 and R2 are equal.
4. Map paired-end FASTQ reads to the sacCer2 genome using Bowtie2.
5. Convert the SAM output to sorted and indexed BAM files.
6. Compute statistics: total, mapped, and unmapped read pairs.
7. Save individual mapping stats and append results to a summary CSV.

---

## 🧪 Tools Required

- `samtools`
- `bedtools2`
- `bowtie2`

These tools must be installed or loaded via environment modules.

---

## 🧬 Input & Output

- **Input**:
  - BAM files containing reads mapped to SK1 chromosomes:
    ```
    /uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/6_mapped_SK1_reads/*.bam
    ```
  - sacCer2 reference index:
    ```
    /uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/8_sacCer2_Referance_genome/bowtie2_index_sacCer2/sacCer2_combined
    ```

- **Output**:
  - Mapped BAM files, stats, and logs in:
    ```
    /uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/9_sk1_reads_mapped_to_sacCer2/
    ```

  - Summary CSV:
    ```
    mapping_summary.csv
    ```

---

## 🖥️ SLURM Script Explained

```bash
#!/bin/bash
#SBATCH --job-name=sk1_mapped_sacCer2         # SLURM job name
#SBATCH --output=sk1_mapped_sacCer2.log       # File for standard output
#SBATCH --error=sk1_mapped_sacCer2.err        # File for standard error
#SBATCH --ntasks=1                            # Number of tasks
#SBATCH --cpus-per-task=8                     # Number of CPU cores
#SBATCH --mem=32G                             # Memory limit
#SBATCH --time=1-00:00:00                     # Time limit (1 day)
#SBATCH --mail-type=BEGIN,END,FAIL            # Send emails on job events
#SBATCH --mail-user=a.rajalingam.23@abdn.ac.uk # Email for notifications

set -euo pipefail                             # Safe script execution

# Load required software
module load samtools
module load bedtools2
module load bowtie2

# Define paths
bam_dir="/.../6_mapped_SK1_reads"             # Input BAM files
output_dir="/.../9_sk1_reads_mapped_to_sacCer2"
ref_index="/.../8_sacCer2_Referance_genome/bowtie2_index_sacCer2/sacCer2_combined"

mkdir -p "$output_dir"                        # Create output directory if missing

summary_csv="${output_dir}/mapping_summary.csv"
echo "Sample,Input Reads (Pairs),Mapped Reads (Pairs),Unmapped Reads (Pairs)" > "$summary_csv"

samples=(...)                                 # List of sample names

for sample in "${samples[@]}"; do
    bam_input="${bam_dir}/${sample}_mapped_SK1.bam"
    sorted_bam="${output_dir}/${sample}_sorted_by_name.bam"
    r1_fastq="${output_dir}/${sample}_R1.fastq"
    r2_fastq="${output_dir}/${sample}_R2.fastq"

    sam_output="${output_dir}/${sample}_mapped_to_sacCer2.sam"
    bam_output="${output_dir}/${sample}_mapped_to_sacCer2.bam"
    sorted_bam_final="${output_dir}/${sample}_mapped_to_sacCer2_sorted.bam"
    mapping_log="${output_dir}/${sample}_mapping.log"
    stats_file="${output_dir}/${sample}_mapping_stats.txt"

    # Sort BAM by read name for paired FASTQ conversion
    samtools sort -n -o "$sorted_bam" "$bam_input"

    # Convert to paired FASTQ
    bedtools bamtofastq -i "$sorted_bam" -fq "$r1_fastq" -fq2 "$r2_fastq"

    # Count lines in FASTQ to validate paired reads
    r1_count=$(wc -l < "$r1_fastq")
    r2_count=$(wc -l < "$r2_fastq")
    if [[ "$r1_count" -ne "$r2_count" ]]; then
        echo "Error: Mismatched read counts for $sample"
        continue
    fi

    # Bowtie2 mapping
    bowtie2 -x "$ref_index" -1 "$r1_fastq" -2 "$r2_fastq" -S "$sam_output"             -N 1 --no-discordant --no-mixed -X 1000 --threads $SLURM_CPUS_PER_TASK             2> "$mapping_log"

    # Convert and sort mapped reads
    samtools view -bS "$sam_output" > "$bam_output"
    samtools sort -o "$sorted_bam_final" "$bam_output"
    samtools index "$sorted_bam_final"

    # Clean up and compress FASTQ
    rm "$sam_output" "$bam_output" "$sorted_bam"
    gzip "$r1_fastq" "$r2_fastq"

    # Calculate stats
    total_reads=$((r1_count / 4))
    mapped_reads=$(samtools view -c -f 2 "$sorted_bam_final")
    mapped_read_pairs=$((mapped_reads / 2))
    unmapped_reads=$((total_reads - mapped_read_pairs))

    {
        echo "Sample: $sample"
        echo "Input Reads (Pairs): $total_reads"
        echo "Mapped Reads (Pairs, Proper Orientation): $mapped_read_pairs"
        echo "Unmapped Reads (Pairs): $unmapped_reads"
    } > "$stats_file"

    echo "$sample,$total_reads,$mapped_read_pairs,$unmapped_reads" >> "$summary_csv"
done

echo "$(date): All SK1 reads successfully mapped to sacCer2."
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
