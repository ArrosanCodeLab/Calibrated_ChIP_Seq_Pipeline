# Step 4: Build Bowtie2 Index from Combined Reference Genome

This is the fourth step of the **Calibrated ChIP-seq Analysis Pipeline**, where we build a Bowtie2 index for a combined reference genome (`SK1` + `Smik`). This step ensures that your trimmed reads can be properly mapped to a merged genome sequence using Bowtie2.

---

## 🛠️ What Does This Step Do?

- Loads Bowtie2 (if necessary)
- Creates a directory for the output Bowtie2 index
- Builds the index from a combined `.fa` file using `bowtie2-build`

---

## 🚀 Prerequisites

Before running this script, make sure the following tools are available:

- `bowtie2` – the aligner and indexing tool

### 🔧 Module Loading (if on HPC)

```bash
module load bowtie2
```

> You may skip this if `bowtie2` is installed on your PATH.

---

## 📂 Folder & File Setup

Before running this step, make sure your working directory has the following structure:

```
/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/
├── 4_reference_genomes/
│   ├── SK1/
│   ├── Smik/
│   ├── genome_seq_SK1_and_Smik.fa
```

> ✅ **Note**: If you cloned the GitHub repo correctly, this structure will already be in place.

The combined genome FASTA file is:
```
/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/4_reference_genomes/genome_seq_SK1_and_Smik.fa
```

The Bowtie2 index will be saved to:
```
/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/4_reference_genomes/combined_genome_index/
```

---

## 📜 Full Bowtie2 Indexing Script

```bash
#!/bin/zsh

# ===============================================================
# 🧬 Script: build_bowtie2_index.zsh
# 🔧 Purpose: To build a Bowtie2 index from the combined SK1 + Smik genome
# 📦 Required: Bowtie2 installed and accessible via command line
# 🔄 Input: Combined genome FASTA file (SK1 + Smik)
# 🧾 Output: Bowtie2 index files (.bt2) in a dedicated directory
# ===============================================================

# ======================= IMPORTANT ==============================
# 🔁 Before using this script:
# 1. Ensure that you have downloaded the folder `3_reference_genomes`
#    from the GitHub repository and placed it in:
#    /uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq
# 2. That folder must contain:
#    - genome_seq_SK1_and_Smik.fa
#    - SK1/ and Smik/ folders with original references (optional)
# ==============================================================

# ===> Path to the combined genome FASTA file <===
genome_fasta="/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/4_reference_genomes/genome_seq_SK1_and_Smik.fa"

# ===> Output directory where Bowtie2 index will be stored <===
output_index="/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/4_reference_genomes/combined_genome_index"

# ===> Create the output directory if it doesn’t already exist <===
mkdir -p "${output_index}"

# ===> Run Bowtie2 to build index <===
echo "🔨 Building Bowtie2 index for ${genome_fasta}..."
bowtie2-build -f "${genome_fasta}" "${output_index}/genome_seq_SK1_and_Smik"

# ===> Check for successful completion and notify user <===
if [ $? -eq 0 ]; then
    echo "✅ Bowtie2 index building completed successfully."
else
    echo "❌ Error: Bowtie2 index building failed." >&2
fi
```

---

## ✅ Checklist Before Running

| Requirement                                                          | Ready? |
|----------------------------------------------------------------------|--------|
| Combined genome FASTA file in correct path                          | ⬜️     |
| Correct folder structure under `3_reference_genomes/`               | ⬜️     |
| Enough storage space (approx. 4–5x the genome file size)            | ⬜️     |
| Script permission (`chmod +x build_bowtie2_index.zsh`)              | ⬜️     |
| `bowtie2` installed and available on path or via module             | ⬜️     |

---

## 🔁 Next Step

➡️ Proceed to: **Step 5 – Mapping Trimmed Reads to Combined Genome using Bowtie2**

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
