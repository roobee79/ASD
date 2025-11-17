# ASD

### 📦 Environment Setup
```bash
conda env create -f asd_gat.yml
conda activate asd_gat
```

## 0. Resource Preparation

Before running the pipeline, download the required Hail resource tables (.ht)
and store them under the following directory:

`Resources/`

1️⃣ Install gsutil
`conda install -c conda-forge gsutil`

2️⃣ Download Required Resources
(1) gnomAD v3.1.1 sites (for variant-level AF filtering)

`gsutil -m cp -r gs://gcp-public-data--gnomad/release/3.1.1/ht/genomes/gnomad.genomes.v3.1.1.sites.ht ./`

(2) gnomAD v2.1.1 LoF metrics (for pLI and LOEUF)

`gsutil -m cp -r gs://gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_transcript.ht ./`

(3) MPC38.ht

You can download the original MPC data from the Broad Institute FTP:
```bash
https://ftp.broadinstitute.org/pub/ExAC_release/release1/regional_missense_constraint/
```

The provided data correspond to hg19/GRCh37 coordinates.
To use it with this pipeline (GRCh38 reference), perform liftover to hg38.

```
🧭 Conversion Example (Hail)
import hail as hl

# Import MPC text file (from FTP)
mpc = hl.import_table("MPC.txt.bgz", impute=True)

# Add locus field (hg19)
mpc = mpc.annotate(
    locus = hl.locus(mpc.chrom, hl.int(mpc.pos), reference_genome='GRCh37')
)

# Liftover to GRCh38
mpc = hl.liftover(mpc, "gs://hail-common/references/grch37_to_grch38.over.chain.gz", reference_genome='GRCh38')

# Key by locus and write Hail Table
mpc = mpc.key_by(mpc.locus)
mpc.write("MPC38.ht", overwrite=True)
```



(4) LCR-hs38.bed

```bash
wget https://github.com/lh3/varcmp/raw/master/scripts/LCR-hs38.bed.gz
gunzip LCR-hs38.bed.gz
```

Then, organize the folder structure as follows:
```bash
Resources/
├── gnomad.genomes.v3.1.1.sites.ht/
├── gnomad.v2.1.1.lof_metrics.by_transcript.ht/
├── MPC38.ht
└── LCR-hs38.bed/                        
```

-----------------------------------------------------------------------------------------------------------

Step 1,2 are modularized in `vcf_process_func.py`

## 1. `filter_oriv.py` - ORIV Filtering and Annotation Pipeline
### 📘 Overview
This code provides a reproducible **Hail-based pipeline** for filtering, annotating, and classifying  
**One-sided Rare Inherited Variants (ORIV)** in human cohort sequencing data.

### 🚀 How to Use
1️⃣ Edit configuration in filter_oriv.py

Important:
You must edit the Configuration and File Path Definitions sections before running the script.

Open the script and set up:
```bash
# Configuration section
date = "run_date"                        # Run date or tag
project = "my_project"                 # Project name
i_dir = "/absolute/path/to/main_dir"   # Main working directory

# File path definitions
vcf_path = "/path/to/raw.vcf.gz"
ped_path = "/path/to/pedigree.ped"
sample_info_path = "/path/to/sample_info.tsv"
```




2️⃣ Sample Information File (sample_info.tsv)

This file provides sample-level metadata used during the Sample_QC() step
to annotate each sample with family and role information.

| Column | Description                                          |
| ------ | ---------------------------------------------------- |
| `s`    | Sample ID (unique per individual)                    |
| `ROLE` | Relationship role within family (`d`, `m`, `p`, `s`) |
| `fam`  | Family ID linking related samples                    |

ROLE legend:

d — Father (dad)

m — Mother

p — Proband (case / affected child)

s — Sibling (unaffected child)

📄 Example Format
```bash
s       ROLE    fam
iid12345       m       14109
iid12346       d       14109
iid12347       p       14109
iid12348       s       14109
```


3️⃣ Run the pipeline

From the repository root:
```bash
python scripts/filter_oriv.py
```

This will sequentially perform:
1. **VCF import and sample QC**
2. **Variant QC (multi-allelic, LCR, call rate, HWE)**
3. **VEP annotation** with LOFTEE, CADD, and MPC
4. **Rare heterozygous filtering**
5. **One-sided rare inherited variant (ORIV) detection**
6. **Export for deep learning / downstream statistical models**

   **Processed Hail MatrixTables (.mt, .ht)**

   **Summary files (.tsv.bgz) for analysis**



All checkpoints and outputs will be saved under:
```bash
{i_dir}/{project}/Inputs/
{i_dir}/{project}/Outputs/
```


#### 📤 Output Files
| Step           | Output                                 | Description                      |
| -------------- | -------------------------------------- | -------------------------------- |
| Sample QC      | `<project>_<date>_after_sampleQC.mt`   | After sample annotation          |
| Variant QC     | `<project>_<date>_after_variantQC.mt`  | After variant-level filtering    |
| VEP            | `<project>_<date>_after_vep.mt`        | After annotation                 |
| Rare Filtering | `<project>_<date>_hqrarehet.mt`        | High-quality rare hets           |
| OIV            | `<project>_<date>_oih.ht`              | One-sided inherited variants     |
| ORIV           | `<project>_<date>_orih.ht`             | Rare inherited variants          |
| ORIV Input     | `<project>_orih_input.tsv.bgz`         | Deep learning model input        |
| ORIV Info      | `<project>_orih_<date>.tsv.bgz`        | Full variant annotation table    |






## 2. `make_akoriv.py` — Cross-Cohort ORIV Integration Pipeline
### 📘 Overview

make_akoriv.py performs cross-cohort inherited variant (ORIV) integration between an internal Korean cohort and an external dataset using the Hail framework.
It identifies shared and unique inherited variants, generates gene-level summaries, and exports model-ready tables for downstream analysis.

This script is designed as the second stage of the variant processing pipeline, following `filter_oriv.py`.

🧩 Input Dependencies

This script requires outputs from filter_oriv.py:

| File Type                | Example Path                                              | Description                              |
| ------------------------ | --------------------------------------------------------- | ---------------------------------------- |
| Korean ORIH table        | `/path/to/<project>/Inputs/<project>_<date>_orih.ht`      | Annotated rare inherited variants        |
| Korean OIH table         | `/path/to/<project>/Inputs/<project>_<date>_oih.ht`       | Raw inherited variant candidates         |
| ORIH input list          | `/path/to/<project>/Outputs/<project>_orih_input.tsv.bgz` | Gene list used for filtering             |
| External cohort IH table | `/path/to/<project2>/Inputs/<project2>_<date>_oih.ht`     | External dataset inherited variant table |


📤 Output Files
| Output                                        | Path                  | Description                                    |
| --------------------------------------------- | --------------------- | -----------------------------------------------|
| `*_ih_in_<project>orihgene.ht`                | `<project2>/Inputs/`  | IH variants overlapping <project> ORIH genes   |
| `*_fam_oih.ht`                                | `<project2>/Inputs/`  | Family-level inherited variant table           |
| `*_uniq_oriv.ht`                              | `<project2>/Inputs/`  | Unique (non-overlapping) inherited variants    |
| `*_overlap_<project>oriv.ht`                  | `<project2>/Inputs/`  | Shared inherited variants between cohorts      |
| `*_akoriv_in_<project>orihgene.ht`            | `<project2>/Inputs/`  | Combined integrated AKORIV dataset             |
| `*_akoriv_in_<project>orihgene_input.tsv.bgz` | `<project2>/Outputs/` | Model-ready per-sample × gene matrix           |
| `*_akoriv_in_<project>orihgene_tidy.ht`       | `<project2>/Inputs/`  | Cleaned annotated Hail Table                   |
| `*_akoriv_in_<project>orihgene.tsv.bgz`       | `<project2>/Outputs/` | Final summarized variant information           |


🚀 Usage
```bash
python make_akoriv.py
```

Before running, update the following fields inside the script:

```bash
date = "date"
project = 'project_name'
project2 = "external_cohort_project_name"
i_dir = "/path/to/maindir"
sample_info_path = "/path/to/sampleinfo.csv"
```


-----------------------------------------------------------------------------------------------------------


## 3. AutismGAT_stable.py — Downstream GAT Model Training

### 📘 Overview
This code corresponds to **project1**, implementing the downstream **Graph Attention Network (GAT)** training stage  
for modeling **rare inherited variant (RIV)** profiles in ASD cohorts.  

After generating AKORIV tables with `make_akoriv.py`, this step performs **graph-based learning** of gene–gene relationships  
and individual variant burdens to derive **gene-level embeddings** for downstream prediction and biological interpretation.

A pretrained model checkpoint (`AutismGAT.pt`) trained on **project1** is also provided for reproducibility,  
enabling users to reproduce inference and visualization without retraining.

---

### 🔬 Step 1 — Filter Genes

Filter the integrated AKORIV dataset to include only the stable genes used for model training:

```bash
models/orih_entrz_v2_scgene_L1_INT.csv
```

This file provides the curated list of genes (Entrez IDs) included in the GAT model input.

## 🧩 Dataset Preparation — Input File Formats

Below are anonymized examples of required input files for constructing the **GAT dataset**.  
Each file must be placed under the appropriate `raw/` or `input/` directory before running the preprocessing pipeline.

---

### (1) `st12_gene_ns.csv` — Gene Network Edges  
This file defines **gene–gene connections** (from STRING or BioGRID).  
All gene identifiers should be converted to **Entrez IDs** for consistent graph construction.

| Index | gene_id1 | gene_id2 |
|-------|-----------|----------|
| 0 | 381 | 4074 |
| 1 | 381 | 1595 |
| 2 | 381 | 79090 |
| 3 | 381 | 57414 |

> 🧠 *Tip:* You can build this file from STRING (v12) or BioGRID and map Ensembl IDs to **Entrez IDs** to match the gene-level features used by the model.

---

### (2) `Korean_ASD_WGS_v4_2186.ped` — Family Structure (PED) File  
Standard **PED format** describing family and individual relationships (trio, sibling, etc.).  
Each row represents a single individual.

| fam_id | s | f_id | m_id | sex | pheno |
|---------|---|------|------|------|--------|
| FAM001 | SAMPLE001 | 0 | 0 | 1 | 1 |
| FAM001 | SAMPLE002 | 0 | 0 | 2 | 1 |
| FAM001 | SAMPLE003 | SAMPLE001 | SAMPLE002 | 1 | 2 |
| FAM002 | SAMPLE004 | 0 | 0 | 1 | 1 |
| FAM002 | SAMPLE005 | 0 | 0 | 2 | 1 |
| FAM002 | SAMPLE006 | SAMPLE004 | SAMPLE005 | 2 | 2 |

> 🔹 `sex`: 1 = male, 2 = female  
> 🔹 `pheno`: 1 = unaffected (control), 2 = affected (ASD case)

---

### (3) `ko_info_ac2.csv` — Sample Annotation Metadata  
Contains per-sample metadata including ASD affection status and optional clustering scores.  
Even if `module` and `css` values are unknown, **the columns must still be present** (use `NaN` or leave empty).


| Index | s | AC_check | module | css |
|--------|---|-----------|--------|-----|
| 0 | SAMPLE001 | 1 | 1 | 9 |
| 1 | SAMPLE002 | 1 | 1 | 8 |
| 2 | SAMPLE003 | 1 | 1 | 6 |
| 3 | SAMPLE004 | 1 | 2 | 10 |
| 4 | SAMPLE005 | 1 | NaN | NaN |
| 5 | SAMPLE006 | 1 | 2 | 6 |


> ⚙️ **Column definitions:**  
> - `s`: Sample ID (unique individual ID)  
> - `AC_check`: ASD affection status (`1` = control, `2` = ASD case)  
> - `module`: Optional numeric label for clustering or network module (must exist, can be `NaN`)  
> - `css`: Optional continuous score (must exist, can be `NaN`)  

> ⚠️ *Note:* The CSV must contain **exactly 5 columns** (`Index`, `s`, `AC_check`, `module`, `css`) — even if `module` and `css` values are missing.  
> The parser expects a fixed schema for proper alignment during data loading.





🚀 Step 2 — Run the GAT Model

Execute the stable GAT model using:

```bash
python models/AutismGAT_stable.py
```

This script loads the filtered AKORIV dataset and trains the Graph Attention Network (GAT)
to learn gene-level embeddings for rare inherited variant profiles across cohorts.


If you wish to use the pretrained model instead of retraining:

```bash
import torch
from models import AutismGAT_stable

model = AutismGAT_stable()
model.load_state_dict(torch.load("models/AutismGAT.pt"))
model.eval()
```

🧩 The provided AutismGAT.pt corresponds to the best-performing model trained on project1,
allowing direct inference, visualization (e.g., attention maps), and downstream analysis without retraining.

-----------------------------------------------------------------------------------------------------------
### 🔐 Data Availability (MSSNG Controlled Access)

The following datasets are currently undergoing MSSNG’s formal release process:

- Step 1 final outputs

- Step 2 final outputs

- Step 3 input data and model weights

These materials will become accessible through the MSSNG controlled-access portal upon completion of the release process.
