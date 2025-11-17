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

(3) LCR-hs38.bed

```bash
wget https://github.com/lh3/varcmp/raw/master/scripts/LCR-hs38.bed.gz
gunzip LCR-hs38.bed.gz
```

Then, organize the folder structure as follows:
```bash
Resources/
├── gnomad.genomes.v3.1.1.sites.ht/
├── gnomad.v2.1.1.lof_metrics.by_transcript.ht/
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
You must edit the Configuration and File Path Definitions sections before running the script.




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



All checkpoints and outputs will be saved under:
```bash
{i_dir}/{project}/Inputs/
{i_dir}/{project}/Outputs/
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
After generating AKORIV tables with make_akoriv.py,
this step performs graph-based modeling of rare inherited variant profiles using a Graph Attention Network (GAT).
The model learns gene-level embeddings from the filtered AKORIV dataset for downstream prediction or interpretation.

🔬 Step 1 — Filter Genes

Filter the integrated AKORIV dataset to include only the stable genes used for model training:

```bash
models/orih_entrz_v2_scgene_L1_INT.csv
```

This file provides the curated list of genes (Entrez IDs) included in the GAT model input.

🚀 Step 2 — Run the GAT Model

Execute the stable GAT model using:

```bash
python models/AutismGAT_stable.py
```

This script loads the filtered AKORIV dataset and trains the Graph Attention Network (GAT)
to learn gene-level embeddings for rare inherited variant profiles across cohorts.


-----------------------------------------------------------------------------------------------------------
### 🔐 Data Availability (MSSNG Controlled Access)

The following datasets are currently undergoing MSSNG’s formal release process:

- Step 1 final outputs

- Step 2 final outputs

- Step 3 input data and model weights

These materials will become accessible through the MSSNG controlled-access portal upon completion of the release process.
