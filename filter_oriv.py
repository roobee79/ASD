"""
# filter_oriv.py
Main ORIV (One-sided Rare Inherited Variant) Filtering and Annotation Pipeline
Author: ERCSB (2025)
"""

import os
import sys
import subprocess
import glob
import hail as hl
import vcf_process_func as VPF


# ============================================================
# 0. Configuration
# ============================================================

# File naming information
date = "date"
project = "project_name"
i_dir = "/path/to/maindir"

# Temporary directory
temp_dir = os.path.join(i_dir, "tmp")
os.makedirs(temp_dir, exist_ok=True)

# Spark configuration
config = {
    "spark.driver.memory": "200g",
    "spark.local.dir": temp_dir,
}

# Initialize Hail
hl.init(
    spark_conf=config,
    master="local[*]",
    tmp_dir=temp_dir,
    local_tmpdir=temp_dir,
    log=f"{i_dir}/{project}/log/{date}.log",
    default_reference="GRCh38",
)

# ============================================================
# 1. File Path Definitions
# ============================================================

vcf_path = "/path/to/rawvcf.vcf"
mt_dir = os.path.join(i_dir, project, "Inputs", "raw_vcf.mt")
ped_path = "/path/to/pedfile.ped"               # PED file
sample_info_path = "/path/to/sampleinfo.csv"    # Sample info file

bed_path = os.path.join(i_dir, "Resources", "LCR-hs38.bed")
gnomAD_path = os.path.join(i_dir, "Resources", "gnomad.genomes.v3.1.1.sites.ht")
# AF_path = os.path.join(i_dir, "Resources", "Public_DB_outerjoin_AF_table_220719.ht")
pLI_score_path = os.path.join(i_dir, "Resources", "gnomad.v2.1.1.lof_metrics.by_transcript.ht")
MPC_path = "/path/to/MPC38.ht"


# ============================================================
# 2. Import and Preprocess
# ============================================================

print("### Step 0. Import and check dataset ###")
mt = hl.import_vcf(vcf_path, reference_genome="GRCh38", n_partitions=300)
mt.write(mt_dir)

# ============================================================
# 3. Sample Filtering and QC
# ============================================================

mt = hl.read_matrix_table(mt_dir)
print("Before filtering samples, MatrixTable count:", mt.count())

# Filter samples based on PED file
mt = VPF.Filter_Sample(ped_path, mt)
mt_dir = os.path.join(i_dir, project, "Inputs", f"{project}_{date}.mt")
print("Checkpoint:", mt_dir)
mt = mt.checkpoint(mt_dir, overwrite=True)

# Sample QC
mt = VPF.Sample_QC(ped_path, sample_info_path, mt)
sampleqc_dir = os.path.join(i_dir, project, "Inputs", f"{project}_{date}_after_sampleQC.mt")
print("Checkpoint:", sampleqc_dir)
mt = mt.checkpoint(sampleqc_dir, overwrite=True)

# ============================================================
# 4. Entry Field Adjustment
# ============================================================

# Ensure AD array matches allele length
mt = mt.annotate_entries(
    DP=hl.sum(mt.AD),
    AD=hl.if_else(
        hl.is_defined(mt.AD),
        hl.if_else(
            hl.len(mt.AD) == hl.len(mt.alleles),
            mt.AD,
            hl.flatten(
                [mt.AD, hl.range(hl.len(mt.alleles) - hl.len(mt.AD)).map(lambda _: 0)]
            ),
        ),
        hl.range(hl.len(mt.alleles)).map(lambda _: 0),
    ),
)

# ============================================================
# 5. Variant QC
# ============================================================

mt = VPF.Variant_QC(bed_path, mt)
variantqc_dir = os.path.join(i_dir, project, "Inputs", f"{project}_{date}_after_variantQC.mt")
print("Checkpoint:", variantqc_dir)
mt = mt.checkpoint(variantqc_dir, overwrite=True)

# Genotype QC
mt = VPF.Genotype_QC(mt)
write_dir = os.path.join(i_dir, project, "Inputs", f"{project}_{date}_after_QC.mt")
print("Checkpoint:", write_dir)
mt = mt.checkpoint(write_dir, overwrite=True)

# ============================================================
# 6. Variant Annotation and Classification
# ============================================================

# VEP Annotation
mt = VPF.Run_vep(i_dir, mt)
after_vep_dir = os.path.join(i_dir, project, "Inputs", f"{project}_{date}_after_vep.mt")
print("Checkpoint:", after_vep_dir)
mt = mt.checkpoint(after_vep_dir, overwrite=True)

# # Classification
# mt = hl.read_matrix_table(after_vep_dir)
# mt = VPF.Classify_variants(MPC_path, gnomAD_path, mt)
# after_classified_dir = os.path.join(i_dir, project, "Inputs", f"{project}_{date}_after_classified.mt")
# print("Checkpoint:", after_classified_dir)
# mt = mt.checkpoint(after_classified_dir, overwrite=True)

# ============================================================
# 7. Variant Filtering Workflow
# ============================================================

# Additional rare variant QC
mt = VPF.Variant_QC2_rare(mt)

# Compute QD
mt = mt.annotate_rows(
    info=mt.info.annotate(
        QD=hl.if_else(
            hl.agg.sum(mt.DP) > 0,
            mt.qual / hl.agg.sum(mt.DP),
            hl.missing(hl.tfloat64),
        )
    )
)

# Filter high-quality rare heterozygous variants
mt = VPF.Filter_HQ_rare_het(mt)

# Save checkpoint
mt_path = os.path.join(i_dir, project, "Inputs", f"{project}_{date}_hqrarehet.mt")
print("Checkpoint:", mt_path)
mt = mt.checkpoint(mt_path, overwrite=True)

# ============================================================
# 8. ORIV Table Generation
# ============================================================

# Generate one-sided inherited variants (OIH)
tb = VPF.MAKE_OIH(mt)
tb_path = os.path.join(i_dir, project, "Inputs", f"{project}_{date}_oih.ht")
print("Checkpoint:", tb_path)
tb = tb.checkpoint(tb_path, overwrite=True)

# Find rare inherited variants (RIH)
tb = VPF.Find_RIH(tb, n_rih_fam=1)
tb_orih_path = os.path.join(i_dir, project, "Inputs", f"{project}_{date}_orih.ht")
print("Checkpoint:", tb_orih_path)
tb = tb.checkpoint(tb_orih_path, overwrite=True)

# ============================================================
# 9. Deep Learning Input Conversion
# ============================================================

print("### ORIV Input Generation ###")
tb_inp = VPF.convert_input(tb)
tb_dir = os.path.join(i_dir, project, "Outputs", f"{project}_orih_input.tsv.bgz")
print("Export:", tb_dir)
tb_inp.export(tb_dir)
print("@@ Finished saving ORIV input file.\n")

# ============================================================
# 10. Detailed ORIV Variant Information
# ============================================================

# tb = VPF.Categorize_allele_frequency(AF_path, tb)
tb = VPF.Annotate_pLI(pLI_score_path, tb)
tb_tidy = VPF.Tidy_up(tb)

tb_tidy_dir = os.path.join(i_dir, project, "Inputs", f"{project}_{date}_orih_tidy.ht")
print("Checkpoint:", tb_tidy_dir)
tb_tidy = tb_tidy.checkpoint(tb_tidy_dir)

print("### Detailed ORIV Info Export ###")
tb_detailed_info_dir = os.path.join(i_dir, project, "Outputs", f"{project}_orih_{date}.tsv.bgz")
print("Export:", tb_detailed_info_dir)
tb_tidy.export(tb_detailed_info_dir)
print("@@ Finished saving detailed ORIV info file.\n")
