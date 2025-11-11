"""
# make_akoriv.py

SSC AKCORIV Variant Filtering and Integration Pipeline

Description:
    This script performs cross-cohort rare inherited variant (ORIV) filtering 
    between external cohort and original cohorts using the Hail framework.
    It identifies shared and unique inherited variants, generates 
    gene-level summaries, and produces model-ready input tables.

Author: ERCSB
Date: 2025-02-04
"""

import os
import sys
import glob
import subprocess
import datetime as dt
import hail as hl
import pandas as pd
import vcf_process_func as VPF 


# ────────────────────────────────────────────────────────────────
# Configuration
# ────────────────────────────────────────────────────────────────
date = "date"
project = 'project_name'
project2 = "external_cohort_project_name"
i_dir = "/path/to/maindir"
job_name = "akoriv"
temp_dir = os.path.join(i_dir, "tmp")

config = {
    "spark.driver.memory": "500g",
    "spark.local.dir": temp_dir,
}

# ────────────────────────────────────────────────────────────────
# Initialize Hail
# ────────────────────────────────────────────────────────────────
os.makedirs(f"{i_dir}/{project2}/log", exist_ok=True)
hl.init(
    spark_conf=config,
    master="local[*]",
    tmp_dir=temp_dir,
    local_tmpdir=temp_dir,
    log=f"{i_dir}/{project2}/log/{job_name}.log",
    default_reference="GRCh38",
)

# ────────────────────────────────────────────────────────────────
# Define resource paths
# ────────────────────────────────────────────────────────────────
gnomAD_path = os.path.join(i_dir, "Resources", "gnomad.genomes.v3.1.1.sites.ht")
AF_path = os.path.join(i_dir, "Resources", "Public_DB_outerjoin_AF_table_220719.ht")
pLI_score_path = os.path.join(i_dir, "Resources", "gnomad.v2.1.1.lof_metrics.by_transcript.ht")

# ────────────────────────────────────────────────────────────────
# Load clinical sample information
# ────────────────────────────────────────────────────────────────
sample_info_path = '/path/to/sampleinfo.csv'
sample_info = hl.import_table(sample_info_path, key="s", delimiter=",")






# ────────────────────────────────────────────────────────────────
# Load reference Korean ORIH data
# ────────────────────────────────────────────────────────────────
kor_oih = hl.read_table(os.path.join(i_dir, project, "Inputs", f"{project}_{date}_oih.ht"))
kor_oih = kor_oih.annotate(variant=hl.variant_str(kor_oih.locus, kor_oih.alleles)).key_by("variant")

kor_orih = hl.read_table(os.path.join(i_dir, project, "Inputs", f"{project}_{date}_orih.ht"))
kor_orih = kor_orih.annotate(variant=hl.variant_str(kor_orih.locus, kor_orih.alleles)).key_by("variant")


# ────────────────────────────────────────────────────────────────
# Load SSC inherited variant table
# ────────────────────────────────────────────────────────────────
ih_path = os.path.join(i_dir, project2, "Inputs", f"{project2}_{date}_oih.ht")
ih = hl.read_table(ih_path)
ih = ih.annotate(variant = hl.variant_str(ih.locus, ih.alleles))

# ────────────────────────────────────────────────────────────────
# Load Korean ORIH input gene list
# ────────────────────────────────────────────────────────────────
kor_orih_input_path = os.path.join(i_dir, project, "Outputs", f"{project}_orih_input.tsv")
kor_orih_input = pd.read_csv(kor_orih_input_path, sep="\t", usecols=["gene_id"]).drop_duplicates().dropna()
kor_orih_gene = hl.Table.from_pandas(kor_orih_input, key="gene_id")

print(f"# of gene in ih: {len(ih.aggregate(hl.agg.counter(ih.gene_id)))}")
print(f"ih count: {ih.count()}")

# ────────────────────────────────────────────────────────────────
# Filter SSC IH variants overlapping Korean ORIH genes
# ────────────────────────────────────────────────────────────────
ih_in_kororihgene = ih.filter(hl.is_defined(kor_orih_gene[ih.gene_id]))
ih_in_kororihgene_dir = os.path.join(
    i_dir, project2, "Inputs", f"{project2}_ih_in_{project}orihgene.ht"
)
ih_in_kororihgene = ih_in_kororihgene.checkpoint(ih_in_kororihgene_dir, overwrite=True)

print(
    f"# of gene in ih_in_kororihgene: {len(ih_in_kororihgene.aggregate(hl.agg.counter(ih_in_kororihgene.gene_id)))}"
)
print(f"ih_in_kororihgene count: {ih_in_kororihgene.count()}")

# ────────────────────────────────────────────────────────────────
# Create family-level inherited variants
# ────────────────────────────────────────────────────────────────
fam_oih_in_kororihgene = VPF.MAKE_FAM_OIH_tb(ih_in_kororihgene)
fam_oih_in_kororihgene_dir = os.path.join(
    i_dir, project2, "Inputs", f"{project2}_fam_oih.ht"
)
fam_oih_in_kororihgene = fam_oih_in_kororihgene.checkpoint(fam_oih_in_kororihgene_dir, overwrite=True)

print(f"fam_oih_in_kororihgene count: {fam_oih_in_kororihgene.count()}")

# ────────────────────────────────────────────────────────────────
# Separate unique ORIV vs Korean ORIV overlap
# ────────────────────────────────────────────────────────────────
uniq_oriv = fam_oih_in_kororihgene.filter(~hl.is_defined(kor_oih[fam_oih_in_kororihgene.variant]))
uniq_oriv_dir = os.path.join(i_dir, project2, "Inputs", f"{project2}_uniq_oriv.ht")
uniq_oriv = uniq_oriv.checkpoint(uniq_oriv_dir, overwrite=True)

koriv = fam_oih_in_kororihgene.filter(hl.is_defined(kor_orih[fam_oih_in_kororihgene.variant]))
koriv_dir = os.path.join(i_dir, project2, "Inputs", f"{project2}_overlap_{project}oriv.ht")
koriv = koriv.checkpoint(koriv_dir, overwrite=True)

# ────────────────────────────────────────────────────────────────
# Merge both sets → AKORIV
# ────────────────────────────────────────────────────────────────
akcoriv_in_kororihgene = uniq_oriv.union(koriv)
akcoriv_in_kororihgene_dir = os.path.join(
    i_dir, project2, "Inputs", f"{project2}_akoriv_in_{project}orihgene.ht"
)
akcoriv_in_kororihgene = akcoriv_in_kororihgene.checkpoint(akcoriv_in_kororihgene_dir, overwrite=True)

print(f"akcoriv_in_{project}orihgene count: {akcoriv_in_kororihgene.count()}")

# ────────────────────────────────────────────────────────────────
# Generate model input (per-sample, per-gene counts)
# ────────────────────────────────────────────────────────────────
print("### Creating model input ###")
akcoriv_input = VPF.convert_input(akcoriv_in_kororihgene)
akcoriv_input_dir = os.path.join(
    i_dir, project2, "Outputs", f"{project2}_akoriv_in_{project}orihgene_input.tsv.bgz"
)
akcoriv_input.export(akcoriv_input_dir)
print("@@ Finished saving input file.\n")

# ────────────────────────────────────────────────────────────────
# Detailed annotation export
# ────────────────────────────────────────────────────────────────
akcoriv_in_kororihgene = VPF.Categorize_allele_frequency(AF_path, akcoriv_in_kororihgene)
akcoriv_in_kororihgene = VPF.Annotate_pLI(pLI_score_path, akcoriv_in_kororihgene)
akcoriv_in_kororihgene_tidy = VPF.Tidy_up(akcoriv_in_kororihgene)

akcoriv_in_kororihgene_tidy_dir = os.path.join(
    i_dir, project2, "Inputs", f"{project2}_akoriv_in_{project}orihgene_tidy.ht"
)
akcoriv_in_kororihgene_tidy = akcoriv_in_kororihgene_tidy.checkpoint(
    akcoriv_in_kororihgene_tidy_dir, overwrite=True
)

akcoriv_detailed_info_dir = os.path.join(
    i_dir, project2, "Outputs", f"{project2}_akoriv_in_{project}orihgene.tsv.bgz"
)
akcoriv_in_kororihgene_tidy.export(akcoriv_detailed_info_dir)
print("@@ Finished saving detailed info file.\n")

# ────────────────────────────────────────────────────────────────
# Done
# ────────────────────────────────────────────────────────────────
print("=== SSC AKCORIV Variant Integration Completed Successfully ===")
