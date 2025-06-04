# Step 4: Bowtie2 Mapping to Combined SK1 + Smik Genome

This is Step 4 of the **Calibrated ChIP-seq Analysis Pipeline**. It maps trimmed paired-end FASTQ files to a combined reference genome (SK1 + Smik) using Bowtie2 and processes the results with samtools.

---

## 📌 What This Step Does

- Loads `bowtie2` and `samtools` modules
- Maps each sample using `bowtie2`
- Filters out multi-mapping reads and overly distant pairs (>10kb)
- Converts output to BAM format using `samtools`
- Sorts and indexes the final BAM file

---

## 📁 Folder Structure

Make sure your project is structured like this:

```
Calibrated_ChIP_Seq/
├── 2_trimmed_output_q20/                  # Trimmed FASTQ files (*.fastq.gz)
├── 3_reference_genomes/
│   ├── genome_seq_SK1_and_Smik.fa         # Combined genome fasta
│   └── combined_genome_index/             # Bowtie2 index files
└── 4_Combain_genome_mapping_output/       # Output BAM files will be saved here
```

---

## 🛠️ Prerequisites

- Bowtie2 index for `genome_seq_SK1_and_Smik.fa` must be already built using `bowtie2-build`
- `fastq.gz` files should be adapter-trimmed and available in `2_trimmed_output_q20/`
- SLURM system with `bowtie2` and `samtools` modules available

---

## 🧬 Sample Names

Make sure your sample names are correctly defined in the `sample_name` array. These names should match the prefix of your paired-end files. Example:

```
IN817_1_trimmed.fastq.gz and IN817_2_trimmed.fastq.gz
```

---

## 📜 SLURM Script (Full)

```bash
#!/bin/bash
#SBATCH --job-name=bowtie2_mapping             # Job name
#SBATCH --output=bowtie2_mapping_%j.out        # Output file
#SBATCH --error=bowtie2_mapping_%j.err         # Error file
#SBATCH --ntasks=1                             # Number of tasks (single task for Bowtie2)
#SBATCH --cpus-per-task=4                      # CPUs per task
#SBATCH --mem=64G                              # Memory allocation
#SBATCH --time=1-00:00:00                      # Maximum runtime (1 day)
#SBATCH --mail-type=BEGIN,END,FAIL             # Email notifications for job events
#SBATCH --mail-user=a.rajalingam.23@abdn.ac.uk # Email address for notifications

# ======================= IMPORTANT =======================
# Update only the sections marked “==> UPDATE THIS <==”.
# =========================================================

# Load required modules for Bowtie2 and Samtools
module load bowtie2
module load samtools

# ===> UPDATE THIS: Define your sample names <===
# These are the base names for paired-end fastq files (e.g., IN817_1_trimmed.fastq.gz and IN817_2_trimmed.fastq.gz)
sample_name=("IN817" "IN855" "IN819" "IN820" "IP817" "IP855" "IP819" "IP820")

# ===> UPDATE THIS: Path to reference genome index built from genome_seq_SK1_and_Smik.fa <===
index_path="/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/3_reference_genomes/combined_genome_index/genome_seq_SK1_and_Smik"

# ===> UPDATE THIS: Path to trimmed FASTQ files <===
fastq_path="/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/2_trimmed_output_q20"

# ===> UPDATE THIS: Path to output directory <===
output_path="/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/4_Combain_genome_mapping_output"
mkdir -p "${output_path}"  # Create output directory if it doesn’t exist

# ---------------------------------------------------------
# Do NOT modify anything below unless you know what you’re doing.
# ---------------------------------------------------------

# Exit on errors and undefined variables for safe execution
set -euo pipefail

# Loop through each sample to perform Bowtie2 mapping and sorting
for i in "${sample_name[@]}"
do
    echo "🔄 Processing sample: ${i}"

    # Construct paths to input FASTQ files
    read1="${fastq_path}/${i}_1_trimmed.fastq.gz"
    read2="${fastq_path}/${i}_2_trimmed.fastq.gz"
    bam_output="${output_path}/${i}_mapped.bam"

    # Perform mapping with Bowtie2, filter, convert to BAM, sort, and index
    echo "📌 Mapping with Bowtie2 for ${i}..."
    bowtie2 -x "${index_path}"             -1 "${read1}"             -2 "${read2}"             -p 4 -N 0 |         grep -v "XS:" |         awk '($9 <= 10000 && $9 >= -10000) || $1 ~ /^@/' |         samtools view -bS -f 2 -@ 3 |         samtools sort -@ 3 -o "${bam_output}"

    echo "🗂️ Indexing BAM file for ${i}..."
    samtools index "${bam_output}"

    echo "✅ Completed sample: ${i}"
done

echo "🎉 Bowtie2 mapping completed for all samples."
```

---

## ✅ Checklist Before Running

| Requirement                                        | Ready? |
|----------------------------------------------------|--------|
| Bowtie2 index exists in `3_reference_genomes/`     | ⬜️     |
| Input FASTQ files are trimmed and named correctly  | ⬜️     |
| `sample_name` array reflects your actual samples    | ⬜️     |
| SLURM system has `bowtie2` and `samtools` loaded    | ⬜️     |

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

## 📬 Corresponding Author

**Hajime Murakami**  
University of Aberdeen, UK  
📧 hajime.murakami1@abdn.ac.uk
