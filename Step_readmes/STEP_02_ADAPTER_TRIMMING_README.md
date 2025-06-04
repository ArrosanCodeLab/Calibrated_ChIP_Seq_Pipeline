# Step 2: Adapter Trimming using `fastp`

This is the second step of the **Calibrated ChIP-seq Analysis Pipeline**, where raw FASTQ files are cleaned by removing sequencing adapters and low-quality reads.

---

## 🛠️ What This Step Does

- Trims adapters from raw paired-end FASTQ files
- Removes low-quality reads (quality threshold set to Q20)
- Ensures a minimum read length of 20 base pairs
- Outputs high-quality trimmed reads for downstream mapping

---

## 🚀 Prerequisites

Before running this script, make sure the following are available:

- `fastp` installed and available in your environment (or loaded via a module)
- Raw FASTQ files located in `/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/1_Raw_Sequence_Data/`
- Paired-end files follow a naming pattern like `SampleID_1.fastq.gz` and `SampleID_2.fastq.gz`

### 🔧 Installing `fastp` (if not available)

```bash
# Using conda
conda install -c bioconda fastp

# Or using Homebrew (macOS)
brew install fastp
```

---

## 📂 Folder Structure

```
/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/
├── 1_Raw_Sequence_Data/
│   ├── IN817_1.fastq.gz
│   ├── IN817_2.fastq.gz
│   └── ...
└── 2_trimmed_output_q20/
    ├── IN817_1_trimmed.fastq.gz
    ├── IN817_2_trimmed.fastq.gz
    └── ...
```

---

## 🧪 Sample SLURM Script

```bash
#!/bin/bash
#SBATCH --job-name=adapter_trimming             # Job name
#SBATCH --output=adapter_trimming_%j.out        # Output file
#SBATCH --error=adapter_trimming_%j.err         # Error file
#SBATCH --ntasks=1                              # Number of tasks
#SBATCH --cpus-per-task=4                       # Number of CPU cores per task
#SBATCH --mem=16G                               # Memory allocation
#SBATCH --time=04:00:00                         # Max runtime
#SBATCH --mail-type=BEGIN,END,FAIL              # Email notifications for job status
#SBATCH --mail-user=a.rajalingam.23@abdn.ac.uk  # Your email for job notifications

# ======================= IMPORTANT =======================
# Update only the two sections marked “==> UPDATE THIS <==”.
# =========================================================

# Define the output directory for trimmed files
# ===> UPDATE THIS to your desired trimmed output directory <===
output_dir="/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/2_trimmed_output_q20"
mkdir -p "$output_dir"  # Ensure the directory exists

# Define your sample names for paired-end reads
# ===> UPDATE THIS to reflect your actual sample names <===
sample_names=("IN817" "IN855" "IN819" "IN820" "IP817" "IP855" "IP819" "IP820")

# ---------------------------------------------------------
# Do NOT modify anything below unless you know what you're doing.
# ---------------------------------------------------------

# Load the required module for fastp
module load fastp

# Loop over all samples and trim adapters using fastp
# ===> Make sure your FASTQ filenames follow this pattern: ${i}_1.fastq.gz and ${i}_2.fastq.gz
# ===> If your filenames differ, UPDATE the suffix accordingly (e.g., _R1.fastq.gz or .fq.gz)
for i in "${sample_names[@]}"
do
    echo "Processing sample: ${i}"

    fastp -i "/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/1_Raw_Sequence_Data/${i}_1.fastq.gz" \
          -I "/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/1_Raw_Sequence_Data/${i}_2.fastq.gz" \
          -3 \
          -o "$output_dir/${i}_1_trimmed.fastq.gz" \
          -O "$output_dir/${i}_2_trimmed.fastq.gz" \
          -q 20 -l 20

    echo "Completed processing for sample: ${i}"
done

echo "Adapter trimming completed for all samples."
```

---

## ✅ Checklist Before Running

| Requirement                                              | Ready? |
|----------------------------------------------------------|--------|
| Correct folder path in `output_dir`                      | ⬜️     |
| Sample names correctly reflect your dataset              | ⬜️     |
| Input FASTQ file suffixes match (`_1.fastq.gz`)          | ⬜️     |
| `fastp` is installed or available as a module            | ⬜️     |
| Script has execute permission (`chmod +x script.sh`)     | ⬜️     |

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
**Hajime Murakami**  
Murakami Lab  
University of Aberdeen, UK  
📧 hajime.murakami1@abdn.ac.uk

## 🔁 Next Step
➡️ Proceed to: **Step 3 – Read Mapping**