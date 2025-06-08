# Step 6: Extract SK1-Mapped Reads as FASTQ

This step extracts reads that map specifically to **SK1 chromosomes** from previously mapped BAM files. It also generates indexed BAMs and compressed FASTQ files for downstream SK1-specific analysis.

---

## 📌 What This Step Does

- Iterates over mapped BAM files from combined genome alignment
- Filters reads aligned to SK1 chromosomes (`SK1_1` to `SK1_16`)
- Generates a BAM file with only SK1-mapped reads
- Indexes the SK1-specific BAM
- Converts the filtered BAM into a gzipped FASTQ file

---

## 📁 Input & Output Structure

- **Input BAMs** should be in:
  ```
  /uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/5_Combain_genome_mapping_output/
  ```

- **Output** will be generated in:
  ```
  /uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/6_SK1_specific_reads/
  ```

  Each sample will generate:
  - `<sample>_mapped_SK1.bam`
  - `<sample>_mapped_SK1.bam.bai`
  - `<sample>_mapped_SK1.fq.gz`

---

## 🧪 Tools Required

- `samtools` (version >= 1.10 recommended)

---

## 🧬 Sample Names

Samples are defined as:
```
("IN817" "IN855" "IN819" "IN820" "IP817" "IP855" "IP819" "IP820")
```

---

## 📜 SLURM Script (Full)

```bash
#!/bin/bash
#SBATCH --job-name=extract_mapped_SK1
#SBATCH --output=extract_mapped_SK1.out
#SBATCH --error=extract_mapped_SK1.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=a.rajalingam.23@abdn.ac.uk

# Define paths
input_dir="/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/5_Combain_genome_mapping_output"
output_dir="/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/6_SK1_specific_reads"

# Create output directory if not exists
mkdir -p "${output_dir}"

# Array of sample names
sample_names=("IN817" "IN855" "IN819" "IN820" "IP817" "IP855" "IP819" "IP820")

# Extract reads mapped to SK1 chromosomes
for sample in "${sample_names[@]}"; do
    echo "Processing sample: ${sample}"

    input_bam="${input_dir}/${sample}_mapped.bam"
    output_bam="${output_dir}/${sample}_mapped_SK1.bam"
    output_fastq="${output_dir}/${sample}_mapped_SK1.fq.gz"

    # Filter reads mapped to SK1 chromosomes
    echo "Filtering reads mapped to SK1 chromosomes..."
    samtools view -@ 7 -b "${input_bam}" SK1_1 SK1_2 SK1_3 SK1_4 SK1_5 SK1_6 SK1_7 SK1_8 SK1_9 SK1_10 SK1_11 SK1_12 SK1_13 SK1_14 SK1_15 SK1_16 > "${output_bam}"

    # Index the SK1-specific BAM file
    echo "Indexing SK1-specific BAM file..."
    samtools index "${output_bam}"

    # Convert filtered BAM to FASTQ
    echo "Generating FASTQ file..."
    samtools fastq -@ 7 "${output_bam}" | gzip > "${output_fastq}"

    echo "Completed processing for sample: ${sample}"
done

echo "All samples processed successfully!"
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
