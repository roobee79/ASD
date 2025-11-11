#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
# vcf_process_func.py
Hail Variant Processing Functions for Rare Inherited Variant (ORIV) Pipeline
"""

import hail as hl


# -----------------------------
# 1. SAMPLE & VARIANT QC
# -----------------------------
def Sample_QC(ped_path, sample_info_path, mt):
    """Annotate sex and sample role based on PED and sample info tables."""
    print("[Sample QC] Annotating sex and role...")
    fam = hl.import_fam(ped_path)
    mt = mt.annotate_cols(isFemale=fam[mt.s].is_female)
    mt = mt.filter_cols(~hl.is_missing(mt.isFemale))
    sample_info = hl.import_table(sample_info_path, key='s')
    mt = mt.annotate_cols(ROLE=sample_info[mt.s].ROLE, fam=mt.s.split('-')[2][:-1])
    return mt


def Variant_QC(bed_path, mt):
    """Apply basic variant QC: filters, multi-allelic split, and LCR removal."""
    print("[Variant QC] Starting basic QC...")
    print(f"  Before QC: {mt.rows().count()} variants")

    mt = mt.filter_rows(hl.len(mt.filters) == 0)
    mt = hl.split_multi_hts(mt)
    mt = mt.filter_rows(mt.was_split == False)
    lcr_bed = hl.import_bed(bed_path, reference_genome='GRCh38')
    mt = mt.filter_rows(~hl.is_defined(lcr_bed[mt.locus]))
    mt = mt.filter_rows((mt.alleles[0].length() <= 50) & (mt.alleles[1].length() <= 50))

    print(f"  After QC: {mt.rows().count()} variants")
    mt = hl.variant_qc(mt)
    return mt


def Variant_QC2_rare(mt):
    """Apply additional filters for rare variants."""
    print("[Variant QC 2] Applying call rate & HWE filters...")
    mt = mt.filter_rows((mt.variant_qc.call_rate >= 0.1) & (mt.variant_qc.p_value_hwe > 1e-12))
    print(f"  After QC2: {mt.rows().count()} variants")
    return mt


# -----------------------------
# 2. VARIANT ANNOTATION & CLASSIFICATION
# -----------------------------
def Run_vep(i_dir, mt):
    """Run Hail-VEP with config JSON."""
    print("[VEP annotation]")
    vep_dir = f"{i_dir}/Resources/VariantQC_and_Annotation/vep-configuration.condahail.json"
    return hl.vep(mt, vep_dir)


def unrelated_AF(mt):
    """Annotate allele frequency in unrelated samples."""
    unrel_af = hl.variant_qc(mt.filter_cols(mt.s.matches('IBS-ASD-[0-9]+[12]-')))
    mt = mt.annotate_rows(
        unrel_af=unrel_af.index_rows(mt.row_key).variant_qc.AF[1],
        unrel_ac=unrel_af.index_rows(mt.row_key).variant_qc.AC[1],
    )
    return mt


def Classify_variants(MPC_path, gnomAD_path, mt):
    """Classify variants by consequence and MPC/CADD scores."""
    print("[Classify variants]")
    MPC = hl.read_table(MPC_path)
    MPC = MPC.key_by("locus", "alleles", "gene_name")
    mt = mt.annotate_rows(
        gene_symbol=mt.vep.transcript_consequences.gene_symbol[0],
        transcript=mt.vep.transcript_consequences.transcript_id[0],
        gene_id=mt.vep.transcript_consequences.gene_id[0],
    )
    mt = mt.key_rows_by(mt.locus, mt.alleles, mt.gene_symbol)
    mt = mt.annotate_rows(MPC=MPC[mt.locus, mt.alleles, mt.gene_symbol].MPC)
    mt = mt.key_rows_by(mt.locus, mt.alleles)

    gnomad_all = hl.read_table(gnomAD_path)
    non_neuro_gnomad = gnomad_all.annotate(
        AF=gnomad_all.freq[gnomad_all.freq_index_dict['non_neuro-adj']].AF,
        AC=gnomad_all.freq[gnomad_all.freq_index_dict['non_neuro-adj']].AC,
    )
    mt = mt.annotate_rows(
        non_neuro_gnomad_af=non_neuro_gnomad[mt.locus, mt.alleles].AF,
        non_neuro_gnomad_ac=non_neuro_gnomad[mt.locus, mt.alleles].AC,
    )
    mt = unrelated_AF(mt)

    PTV = (
        (hl.array(['frameshift_variant', 'splice_acceptor_variant', 'splice_donor_variant', 'stop_gained'])
         .contains(mt.vep.most_severe_consequence))
        & (mt.vep.transcript_consequences.lof[0] == 'HC')
        & ((mt.vep.transcript_consequences.lof_flags[0] == 'SINGLE_EXON')
           | (hl.is_missing(mt.vep.transcript_consequences.lof_flags[0])))
    )
    damaging_missense = (
        (hl.array(['missense_variant', 'stop_lost', 'start_lost', 'protein_altering_variant'])
         .contains(mt.vep.most_severe_consequence))
        & (((mt.vep.transcript_consequences.polyphen_prediction[0] == 'probably_damaging')
            & (mt.vep.transcript_consequences.sift_prediction[0] == 'deleterious'))
           | (mt.vep.transcript_consequences.cadd_phred[0] > 20)
           | (mt.MPC >= 2))
    )

    mt = mt.annotate_rows(
        VariantType=hl.case(missing_false=True)
        .when(PTV, 'protein truncating')
        .when(damaging_missense, 'damaging missense')
        .or_missing()
    )
    return mt


# -----------------------------
# 3. RARE VARIANT FILTERING
# -----------------------------
def Filter_HQ_rare_het(mt):
    """Filter high-quality rare heterozygous variants."""
    print("[Filter HQ rare heterozygous variants]")
    mt = mt.filter_entries(mt.GT.is_het())
    mt = mt.filter_rows(mt.variant_qc.gq_stats.mean >= 76)
    return mt


def MAKE_OIH(mt):
    """Generate one-sided inherited variant table."""
    mt = mt.filter_rows((mt.non_neuro_gnomad_af <= 0.01) | hl.is_missing(mt.non_neuro_gnomad_af))
    tmt = mt.entries()
    tmt_roleagg = tmt.group_by(tmt.locus, tmt.alleles, tmt.fam).aggregate(
        aggROLE=hl.agg.collect_as_set(tmt.ROLE)
    )
    ih = tmt_roleagg.filter(
        ((tmt_roleagg.aggROLE.contains('p')) | (tmt_roleagg.aggROLE.contains('s')))
        & (tmt_roleagg.aggROLE.contains('f') | tmt_roleagg.aggROLE.contains('m'))
    )
    oih = ih.filter(
        (ih.aggROLE.contains('p') & ~ih.aggROLE.contains('s'))
        | (~ih.aggROLE.contains('p') & ih.aggROLE.contains('s'))
    )
    oih = oih.key_by('locus', 'alleles', 'fam')
    tmt = tmt.key_by('locus', 'alleles', 'fam')
    return tmt.semi_join(oih)


def Find_RIH(tmt_ih, n_rih_fam):
    """Find rare inherited variants (RIH) shared by ≤ n families."""
    tmt_famagg = tmt_ih.group_by(tmt_ih.locus, tmt_ih.alleles).aggregate(
        famset=hl.agg.collect_as_set(tmt_ih.fam)
    )
    ih_fam = tmt_famagg.filter((hl.len(tmt_famagg.famset) <= n_rih_fam) & (hl.len(tmt_famagg.famset) > 0))
    return tmt_ih.semi_join(ih_fam)


def MAKE_FAM_OIH_tb(tmt):
    print('tmt count', tmt.count(),'!')
    tmt_roleagg = tmt.group_by(tmt.locus, tmt.alleles, tmt.fam).aggregate(aggROLE = hl.agg.collect_as_set(tmt.ROLE)) 
    ih = tmt_roleagg.filter(((tmt_roleagg.aggROLE.contains('p'))|(tmt_roleagg.aggROLE.contains('s'))) & (tmt_roleagg.aggROLE.contains('f')|tmt_roleagg.aggROLE.contains('m'))) #p혹은s 중 하나 가지고 있고, 부모 있어야 함, inherited
    ih = ih.filter((ih.aggROLE.contains('p')& ~ih.aggROLE.contains('s')) | (~ih.aggROLE.contains('p')&ih.aggROLE.contains('s'))) # # p만 가지고 있거나 / s만 가지고 있는 경우
    ih = ih.key_by('locus','alleles','fam')
    tmt = tmt.key_by('locus','alleles','fam')
    tmt_ih = tmt.semi_join(ih)  
    tmt = tmt.key_by('locus','alleles')
    print('tmt_ih count', tmt_ih.count(),'!') 
    return tmt_ih

# -----------------------------
# 4. FINAL TABLE GENERATION
# -----------------------------
def Categorize_allele_frequency(AF_path, tb):
    AF_table = hl.read_table(AF_path)
    tb = tb.key_by('locus', 'alleles')
    tb = tb.annotate(AF=AF_table[tb.locus, tb.alleles])
    return tb


def Annotate_pLI(pLI_score_path, tb):
    gnomad_resource = hl.read_table(pLI_score_path)
    tb = tb.key_by("gene_symbol", "transcript")
    tb = tb.annotate(
        pLI=gnomad_resource[tb.gene_symbol, tb.transcript].pLI,
        oe_lof_upper=gnomad_resource[tb.gene_symbol, tb.transcript].oe_lof_upper,
    )
    tb = tb.key_by(tb.locus, tb.alleles, tb.s)
    return tb


def convert_input(tb):
    tb = tb.key_by('locus', 'alleles', 's')
    tb = tb.select(tb.gene_id, tb.variant_class, tb.VariantType)
    tb2 = tb.group_by(tb.s, tb.gene_id)
    numcount = tb2.aggregate(
        numsnp=hl.agg.count_where(tb.variant_class == 'snp'),
        numindel=hl.agg.count_where(tb.variant_class == 'indel'),
        numPTV=hl.agg.count_where(tb.VariantType == 'protein truncating'),
    )
    return numcount
