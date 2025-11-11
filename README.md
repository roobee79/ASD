# ASD

## ORIV Filtering and Annotation Pipeline

This repository provides a reproducible **Hail-based pipeline** for filtering, annotating, and classifying  
**One-sided Rare Inherited Variants (ORIV)** in human cohort sequencing data.

---
### 📦 Environment Setup
```bash
conda env create -f asd_gat.yml
conda activate asd_gat
```
---

### 🧩 Overview

All major steps are modularized in `vcf_process_func.py`,  
and orchestrated by the main driver script `filter_oriv.py`.

### 🚀 How to Use
#### 1️⃣ Edit configuration in filter_oriv.py

Open the script and set up:
```bash
# Configuration section
date = "run_date"                        # Run date or tag
project = "my_project"                 # Project name
i_dir = "/absolute/path/to/main_dir"   # Main working directory

# File path definitions
vcf_path = "/path/to/raw.vcf.gz"
ped_path = "/path/to/pedigree.ped"
sample_info_path = "/path/to/sample_info.csv"
MPC_path = "/path/to/MPC38.ht"
```

⚠️ Important:
You must edit the Configuration and File Path Definitions sections
before running the script.
All checkpoints and outputs will be saved under:
```bash
{i_dir}/{project}/Inputs/
{i_dir}/{project}/Outputs/
```

#### 2️⃣ Run the pipeline

From the repository root:
```bash
python scripts/filter_oriv.py
```

This will sequentially perform:
1. **VCF import and sample QC**
2. **Variant QC (multi-allelic, LCR, call rate, HWE)**
3. **VEP annotation** with LOFTEE, CADD, and MPC
4. **Variant classification** (PTV, damaging missense, etc.)
5. **Rare heterozygous filtering**
6. **One-sided rare inherited variant (ORIV) detection**
7. **Export for deep learning / downstream statistical models**
    Processed Hail MatrixTables (.mt, .ht)
    Summary files (.tsv.bgz) for analysis

#### 📤 Output Files
| Step           | Output                                 | Description                      |
| -------------- | -------------------------------------- | -------------------------------- |
| Sample QC      | `<project>_<date>_after_sampleQC.mt`   | After sample annotation          |
| Variant QC     | `<project>_<date>_after_variantQC.mt`  | After variant-level filtering    |
| VEP            | `<project>_<date>_after_vep.mt`        | After annotation                 |
| Classified     | `<project>_<date>_after_classified.mt` | After consequence classification |
| Rare Filtering | `<project>_<date>_hqrarehet.mt`        | High-quality rare hets           |
| OIH            | `<project>_<date>_oih.ht`              | One-sided inherited variants     |
| ORIH           | `<project>_<date>_orih.ht`             | Rare inherited variants          |
| ORIV Input     | `<project>_orih_input.tsv.bgz`         | Deep learning model input        |
| ORIV Info      | `<project>_orih_<date>.tsv.bgz`        | Full variant annotation table    |
