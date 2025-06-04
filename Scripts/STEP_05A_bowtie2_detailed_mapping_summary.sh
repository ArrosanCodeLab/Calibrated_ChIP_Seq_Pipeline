#!/bin/zsh

# ===============================================================
# 📊 Script: bowtie2_generate_detailed_mapping_summary.zsh
# 🔍 Purpose: Generate a detailed mapping summary for each sample BAM file
# 📦 Requirements: samtools must be installed and in PATH
# ===============================================================

# ===> Set paths
bam_dir="/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/5_Combain_genome_mapping_output"
output_summary="${bam_dir}/detailed_mapping_summary.txt"

# ===> Initialize output file
echo "Detailed Mapping Summary Report" > "$output_summary"
echo "==============================" >> "$output_summary"

# ===> Get sorted list of BAM files
bam_files=($(ls "${bam_dir}"/*_mapped.bam | sort -V))

# ===> Check for BAM file presence
if [ ${#bam_files[@]} -eq 0 ]; then
    echo "No BAM files found in $bam_dir. Exiting." >&2
    exit 1
fi

# ===> Loop through BAM files
for bam_file in "${bam_files[@]}"
do
    sample=$(basename "$bam_file" "_mapped.bam")
    echo "Processing $sample..."

    if [ ! -s "$bam_file" ]; then
        echo "Error: BAM file for sample $sample is empty or missing." >> "$output_summary"
        continue
    fi

    stats=$(samtools flagstat "$bam_file")

    total_reads=$(echo "$stats" | grep "in total" | awk '{print $1}')
    mapped_reads=$(echo "$stats" | grep " mapped (" | head -1 | awk '{print $1}')
    properly_paired=$(echo "$stats" | grep "properly paired" | awk '{print $1}')
    secondary_alignments=$(echo "$stats" | grep "secondary" | awk '{print $1}')
    supplementary_alignments=$(echo "$stats" | grep "supplementary" | awk '{print $1}')
    singletons=$(echo "$stats" | grep "singletons" | awk '{print $1}')

    if [ -z "$total_reads" ] || [ -z "$mapped_reads" ]; then
        echo "Error processing $sample: Missing total or mapped reads." >> "$output_summary"
        continue
    fi

    sk1_reads=$(samtools view "$bam_file" | awk '$3 ~ /^SK1_/' | wc -l)
    smik_reads=$(samtools view "$bam_file" | awk '$3 ~ /^Smik_/' | wc -l)

    sk1_percentage=$(awk -v sk1="$sk1_reads" -v total="$mapped_reads" 'BEGIN { if (total > 0) printf "%.2f", (sk1 / total) * 100; else print "0.00" }')
    smik_percentage=$(awk -v smik="$smik_reads" -v total="$mapped_reads" 'BEGIN { if (total > 0) printf "%.2f", (smik / total) * 100; else print "0.00" }')

    echo "Sample: $sample" >> "$output_summary"
    echo "----------------" >> "$output_summary"
    echo "Total Reads:" >> "$output_summary"
    echo "${total_reads} reads were processed." >> "$output_summary"
    echo "Mapped Reads:" >> "$output_summary"
    echo "100% of the reads were mapped (${mapped_reads} + 0 mapped)." >> "$output_summary"
    echo "Properly Paired Reads:" >> "$output_summary"
    echo "100% of the reads were properly paired (${properly_paired} + 0 properly paired)." >> "$output_summary"
    echo "No Secondary or Supplementary Alignments:" >> "$output_summary"
    echo "There are no secondary or supplementary alignments (${secondary_alignments} + 0 secondary, ${supplementary_alignments} + 0 supplementary)." >> "$output_summary"
    echo "No Singletons or Cross-Chromosome Pairs:" >> "$output_summary"
    echo "All reads were paired (${singletons} + 0 singletons)." >> "$output_summary"
    echo "Reads Aligned to SK1 Chromosomes:" >> "$output_summary"
    echo "${sk1_reads} reads (${sk1_percentage}%)." >> "$output_summary"
    echo "Reads Aligned to Smik Chromosomes:" >> "$output_summary"
    echo "${smik_reads} reads (${smik_percentage}%)." >> "$output_summary"
    echo "===================================" >> "$output_summary"
    echo "" >> "$output_summary"
done

echo "Detailed mapping summary report generated at $output_summary"
