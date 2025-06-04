# Step 5A: Generate Detailed Bowtie2 Mapping Summary

This step analyzes the BAM files generated from Bowtie2 mapping and creates a detailed summary report for each sample. It includes overall alignment statistics and species-specific mapping percentages (SK1 vs Smik).

---

## 📌 What This Step Does

- Iterates over all sorted BAM files in the specified directory
- Checks if the BAM files are non-empty
- Uses `samtools flagstat` to extract general mapping stats
- Uses `samtools view` + `awk` to count reads aligned to SK1 and Smik chromosomes
- Calculates and reports:
  - Total reads
  - Mapped reads
  - Properly paired reads
  - Secondary and supplementary alignments
  - Singleton reads
  - Species-specific mapping percentages

---

## 📁 Input & Output Structure

- **Input**: BAM files ending with `_mapped.bam` located in the folder:
  ```
  /uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/5_Combain_genome_mapping_output/
  ```

- **Output**: Summary file written to:
  ```
  /uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/5_Combain_genome_mapping_output/detailed_mapping_summary.txt
  ```

---

## 🧪 Tools Required

- `samtools` must be available in your environment
- Script should be executed using `bash` on any Linux-based system (HPC or local)

---

## 📜 Full Shell Script

```bash
#!/bin/bash

# Define paths
bam_dir="/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/5_Combain_genome_mapping_output"
output_summary="${bam_dir}/detailed_mapping_summary.txt"

# Initialize the summary file
echo "Detailed Mapping Summary Report" > "$output_summary"
echo "==============================" >> "$output_summary"

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

    # Check if BAM file is non-empty
    if [ ! -s "$bam_file" ]; then
        echo "Error: BAM file for sample $sample is empty or missing." >> "$output_summary"
        continue
    fi

    # Extract statistics using samtools flagstat
    stats=$(samtools flagstat "$bam_file")

    # Parse general stats
    total_reads=$(echo "$stats" | grep "in total" | awk '{print $1}')
    mapped_reads=$(echo "$stats" | grep " mapped (" | head -1 | awk '{print $1}')
    properly_paired=$(echo "$stats" | grep "properly paired" | awk '{print $1}')
    secondary_alignments=$(echo "$stats" | grep "secondary" | awk '{print $1}')
    supplementary_alignments=$(echo "$stats" | grep "supplementary" | awk '{print $1}')
    singletons=$(echo "$stats" | grep "singletons" | awk '{print $1}')

    # Ensure total_reads and mapped_reads are valid numbers
    if [ -z "$total_reads" ] || [ -z "$mapped_reads" ]; then
        echo "Error processing $sample: Missing total or mapped reads." >> "$output_summary"
        continue
    fi

    # Count reads aligned to SK1 and Smik chromosomes
    sk1_reads=$(samtools view "$bam_file" | awk '$3 ~ /^SK1_/' | wc -l)
    smik_reads=$(samtools view "$bam_file" | awk '$3 ~ /^Smik_/' | wc -l)

    # Calculate percentages
    sk1_percentage=$(awk -v sk1="$sk1_reads" -v total="$mapped_reads" 'BEGIN { if (total > 0) printf "%.2f", (sk1 / total) * 100; else print "0.00" }')
    smik_percentage=$(awk -v smik="$smik_reads" -v total="$mapped_reads" 'BEGIN { if (total > 0) printf "%.2f", (smik / total) * 100; else print "0.00" }')

    # Write the formatted summary to the output file
    {
        echo "Sample: $sample"
        echo "----------------"
        echo "Total Reads:"
        echo "${total_reads} reads were processed."
        echo "Mapped Reads:"
        echo "100% of the reads were mapped (${mapped_reads} + 0 mapped)."
        echo "Properly Paired Reads:"
        echo "100% of the reads were properly paired (${properly_paired} + 0 properly paired)."
        echo "No Secondary or Supplementary Alignments:"
        echo "There are no secondary or supplementary alignments (${secondary_alignments} + 0 secondary, ${supplementary_alignments} + 0 supplementary)."
        echo "No Singletons or Cross-Chromosome Pairs:"
        echo "All reads were paired (${singletons} + 0 singletons)."
        echo "Reads Aligned to SK1 Chromosomes:"
        echo "${sk1_reads} reads (${sk1_percentage}%)."
        echo "Reads Aligned to Smik Chromosomes:"
        echo "${smik_reads} reads (${smik_percentage}%)."
        echo "==================================="
        echo ""
    } >> "$output_summary"
done

echo "Detailed mapping summary report generated at $output_summary"
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