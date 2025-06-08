# Step 5B: Generate Bowtie2 Mapping Summary Table

This step extracts summary mapping statistics from each Bowtie2-generated BAM file and compiles them into a tab-delimited table. The summary includes counts of total reads and the number of reads mapped to SK1 and Smik chromosomes, along with their percentages.

---

## 📌 What This Step Does

- Iterates over all BAM files in the specified directory
- Extracts total and mapped reads using `samtools flagstat`
- Uses `samtools view` and `awk` to count reads aligned to:
  - SK1 chromosomes (matching `SK1_`)
  - Smik chromosomes (matching `Smik_`)
- Calculates the percentage of mapped reads for each species
- Appends a row to a summary table for each sample

---

## 📁 Input & Output Structure

- **Input**: BAM files ending with `_mapped.bam` in:
  ```
  /uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/5_Combain_genome_mapping_output/
  ```

- **Output**: Tab-separated summary table written to:
  ```
  /uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/5_Combain_genome_mapping_output/mapping_summary_table.txt
  ```

---

## 🧪 Tools Required

- `samtools` installed and accessible via command line
- `bash` or `zsh` shell environment

---

## 📜 Full Shell Script

```bash
#!/bin/bash

# Define paths
bam_dir="/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/5_Combain_genome_mapping_output"
output_table="${bam_dir}/mapping_summary_table.txt"

# Initialize the output table file
echo -e "Sample\tTotal_Reads\tMapped_to_SK1\tMapped_to_Smik\tSK1_Percentage\tSmik_Percentage" > "$output_table"

# Sort BAM files numerically
bam_files=($(ls "${bam_dir}"/*_mapped.bam | sort -V))

# Check if BAM files exist
if [ ${#bam_files[@]} -eq 0 ]; then
    echo "No BAM files found in $bam_dir. Exiting." >&2
    exit 1
fi

# Loop through each BAM file
for bam_file in "${bam_files[@]}"
do
    sample=$(basename "$bam_file" "_mapped.bam") # Extract sample name
    echo "Processing $sample..."

    # Extract statistics using samtools flagstat
    stats=$(samtools flagstat "$bam_file")

    # Parse general stats
    total_reads=$(echo "$stats" | grep "in total" | awk '{print $1}')
    mapped_reads=$(echo "$stats" | grep " mapped (" | head -1 | awk '{print $1}')

    # Ensure total_reads and mapped_reads are valid numbers
    if [ -z "$total_reads" ] || [ -z "$mapped_reads" ]; then
        echo "Error processing $sample: Missing total or mapped reads." >&2
        continue
    fi

    # Count reads aligned to SK1 and Smik chromosomes
    sk1_reads=$(samtools view "$bam_file" | awk '$3 ~ /^SK1_/' | wc -l)
    smik_reads=$(samtools view "$bam_file" | awk '$3 ~ /^Smik_/' | wc -l)

    # Calculate percentages
    sk1_percentage=$(awk -v sk1="$sk1_reads" -v total="$mapped_reads" 'BEGIN { if (total > 0) printf "%.2f", (sk1 / total) * 100; else print "0.00" }')
    smik_percentage=$(awk -v smik="$smik_reads" -v total="$mapped_reads" 'BEGIN { if (total > 0) printf "%.2f", (smik / total) * 100; else print "0.00" }')

    # Append to the table
    echo -e "${sample}\t${total_reads}\t${sk1_reads}\t${smik_reads}\t${sk1_percentage}\t${smik_percentage}" >> "$output_table"
done

echo "Mapping summary table generated at $output_table"
```

---

## 👨‍🔬 Author

**Arrosan Rajalingam**  
PhD Student – Murakami Lab  
University of Aberdeen, UK  
📧 a.rajalingam.23@abdn.ac.uk

## 🤝 Contributor

**Yusuke Tsuruta**  
Tokyo Metropolitan University, Japan

## 📬 Corresponding Author

**Dr. Hajime Murakami**  
University of Aberdeen, UK  
📧 hajime.murakami1@abdn.ac.uk
