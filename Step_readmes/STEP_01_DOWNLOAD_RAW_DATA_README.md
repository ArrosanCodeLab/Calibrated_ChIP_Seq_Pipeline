# Step 1: Download & Validate Raw Sequence Data

This is the first step of the **Calibrated ChIP-seq Analysis Pipeline**, designed to ensure that your raw sequencing files are correctly downloaded and intact. This script:

1. 📦 Downloads raw FASTQ files from URLs listed in a metadata file  
2. ✅ Validates each download using its MD5 checksum  

---

## 🛠️ What Does This Step Do?

- Creates a directory for raw data if it doesn’t exist  
- Reads a metadata file to retrieve filenames, URLs, and MD5 checksums  
- Downloads `.fastq.gz` files (if not already present)  
- Validates the integrity of each file using `md5sum`  

---

## 🚀 Prerequisites

Before running this script, make sure the following tools are installed:

- `wget` – for downloading files  
- `awk` – for extracting fields from the metadata  
- `md5sum` – for validating checksums  

### 🔧 Installation Commands

```bash
# On Ubuntu/Debian
sudo apt update
sudo apt install wget coreutils

# On macOS with Homebrew
brew install wget coreutils
```

---

## 📂 Folder Structure & Metadata File

Inside your working directory, ensure you have a subdirectory to hold the raw FASTQ downloads:

```bash
mkdir -p /uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/1_Raw_Sequence_Data
```

Then place your metadata file in:

```
/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/1_Raw_Sequence_Data/EN****_22samples_md5sum_DownloadLink.txt
```

> ⚠️ **Important:** If you want to run this pipeline for a different batch of samples, update both:
> 1. `raw_data_dir` (the folder path)  
> 2. `md5sum_file` (the metadata filename)  
>
> For example, if your new folder is `/.../New_Project_Data`, and your metadata file is `New_Metadata.txt`, set:
> ```bash
> raw_data_dir="/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/1_Raw_Sequence_Data/New_Project_Data"
> md5sum_file="${raw_data_dir}/New_Metadata.txt"
> ```

### 📝 Metadata File Format

This file must be **tab-delimited** with at least three columns in this order:

```
Index   Filename            Checksum                             URL
1       Sample_1.fastq.gz   a0c8e3c8120e4c985b98758f8c15e097     https://ftp.url/Sample_1.fastq.gz
2       Sample_2.fastq.gz   93a8d9c762a4d8fc3448f1d456a890bb     https://ftp.url/Sample_2.fastq.gz
...
```

- **Column 1** (`Index`) can be any integer and is ignored by the script (it just skips the header using `NR>1`).  
- **Column 2** (`Filename`) is the name you expect after downloading (e.g., `Sample_1.fastq.gz`; used for checksum validation).  
- **Column 3** (`Checksum`) is the MD5 hash string for the file.  
- **Column 4** (`URL`) is the full HTTP or FTP link to download the `.fastq.gz`.  

> **⚠️ Make sure there is no header**, or if there is, the script skips it via `NR>1` in each `awk` command.

---

## 📜 Line-by-Line Script Breakdown

Below is the exact script. **Do not modify any lines**—instead, update only the variables indicated (e.g., `raw_data_dir`, `md5sum_file`) when running with a new dataset.

```bash
#!/bin/bash
#SBATCH --job-name=download_fastq
#SBATCH --output=download_fastq_%j.out
#SBATCH --error=download_fastq_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=a.rajalingam.23@abdn.ac.uk

# ======================= IMPORTANT =======================
# Update only the two lines marked “==> UPDATE THIS <==”.
# =========================================================

# Define the directory for raw data downloads
# ===> we want to create 1_Raw_Sequence_Data under the main folder <===
raw_data_dir="/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/1_Raw_Sequence_Data"

# Define the path to the metadata file (MD5 + URLs)
# ===> metadata is located one level up, in the main Calibrated_ChIP_Seq folder <===
md5sum_file="/uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/EN****_22samples_md5sum_DownloadLink.txt"

# ---------------------------------------------------------
# Do NOT modify anything below unless you know what you’re doing.
# ---------------------------------------------------------

# Ensure the raw data directory exists (will create 
# /uoa/home/r02ar23/sharedscratch/Calibrated_ChIP_Seq/1_Raw_Sequence_Data)
mkdir -p "$raw_data_dir"

# Change to the raw data directory
cd "$raw_data_dir" || { echo "Cannot change to raw data directory!"; exit 1; }

# Download each file listed in the metadata file
echo "Starting downloads..."
awk 'NR>1 {print $4}' "$md5sum_file" | while read -r url; do
    filename=$(basename "$url")
    if [ -f "$filename" ]; then
        echo "File $filename already exists. Skipping download."
    else
        echo "Downloading $url..."
        wget "$url" || { echo "Failed to download $url"; exit 1; }
    fi
done

# Validate downloads using MD5 checksums
echo "Validating downloads with MD5 checksums..."
awk 'NR>1 {print $3, $1}' "$md5sum_file" | while read -r checksum filename; do
    if md5sum -c <(echo "$checksum $filename"); then
        echo "$filename: VALID"
    else
        echo "$filename: INVALID"
        exit 1
    fi
done

echo "All files downloaded and validated successfully!"
```

---

## ✅ Checklist Before Running

| Requirement                                              | Ready? |
|----------------------------------------------------------|--------|
| Correct folder path in `raw_data_dir`                     | ⬜️     |
| Metadata file located at `${raw_data_dir}/…`               | ⬜️     |
| Metadata columns: Index (or Filename), Filename, Checksum, URL   | ⬜️     |
| Tools installed: `wget`, `awk`, `md5sum`                   | ⬜️     |
| Script has execute permission (`chmod +x download_fastq.sh`)| ⬜️     |
| Internet connection active                                 | ⬜️     |

---

## 🔁 Next Step

➡️ Proceed to: **Step 2 – Adapter Trimming**, where you will remove sequencing adapters from the raw FASTQ files before alignment.

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
➡️ Proceed to: **Step 2 – Adapter Trimming**