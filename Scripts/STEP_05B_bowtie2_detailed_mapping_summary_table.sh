#!/bin/bash

# ================================================================
# 📊 Script: generate_mapping_summary_table.sh
# 🔍 Purpose: Generate a tab-delimited summary table from mapped BAM files
# 📦 Required: samtools must be installed and accessible via the command line
# 🔄 Input: Sorted and indexed BAM files (one per sample)
# 📄 Output: mapping_summary_table.txt with per-sample SK1/Smik alignment stats
# ================================================================

# ===> Define directory containing BAM files <===
bam_dir="/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/5_Combain_genome_mapping_output"

# ===> Define output file for the summary table <===
output_table="${bam_dir}/mapping_summary_table.txt"

# ===> Initialize the output table with a header <===
echo -e "Sample\tTotal_Reads\tMapped_to_SK1\tMapped_to_Smik\tSK1_Percentage\tSmik_Percentage" > "$output_table"

# ===> Get sorted list of BAM files <===
bam_files=($(ls "${bam_dir}"/*_mapped.bam | sort -V))

# ===> Exit if no BAM files are found <===
if [ ${#bam_files[@]} -eq 0 ]; then
    echo "No BAM files found in $bam_dir. Exiting." >&2
    exit 1
fi

# ===> Loop through each BAM file and summarize mapping statistics <===
for bam_file in "${bam_files[@]}"
do
    sample=$(basename "$bam_file" "_mapped.bam")  # Extract sample name
    echo "Processing $sample..."

    # Run samtools flagstat to get general mapping statistics
    stats=$(samtools flagstat "$bam_file")

    # Extract total reads and mapped reads
    total_reads=$(echo "$stats" | grep "in total" | awk '{print $1}')
    mapped_reads=$(echo "$stats" | grep " mapped (" | head -1 | awk '{print $1}')

    # Sanity check for numeric values
    if [ -z "$total_reads" ] || [ -z "$mapped_reads" ]; then
        echo "Error processing $sample: Missing total or mapped reads." >&2
        continue
    fi

    # Count how many reads are aligned to SK1 and Smik
    sk1_reads=$(samtools view "$bam_file" | awk '$3 ~ /^SK1_/' | wc -l)
    smik_reads=$(samtools view "$bam_file" | awk '$3 ~ /^Smik_/' | wc -l)

    # Calculate alignment percentages
    sk1_percentage=$(awk -v sk1="$sk1_reads" -v total="$mapped_reads" 'BEGIN { if (total > 0) printf "%.2f", (sk1 / total) * 100; else print "0.00" }')
    smik_percentage=$(awk -v smik="$smik_reads" -v total="$mapped_reads" 'BEGIN { if (total > 0) printf "%.2f", (smik / total) * 100; else print "0.00" }')

    # Write per-sample stats to the output file
    echo -e "${sample}\t${total_reads}\t${sk1_reads}\t${smik_reads}\t${sk1_percentage}\t${smik_percentage}" >> "$output_table"
done

echo "✅ Mapping summary table generated at $output_table"
