# ASD

## ORIV Filtering and Annotation Pipeline

This repository provides a reproducible **Hail-based pipeline** for filtering, annotating, and classifying  
**One-sided Rare Inherited Variants (ORIV)** in human cohort sequencing data.

---

### 🧩 Overview

The pipeline performs:

1. **VCF import and sample QC**
2. **Variant QC (multi-allelic, LCR, call rate, HWE)**
3. **VEP annotation** with LOFTEE, CADD, and MPC
4. **Variant classification** (PTV, damaging missense, etc.)
5. **Rare heterozygous filtering**
6. **One-sided rare inherited variant (ORIV) detection**
7. **Export for deep learning / downstream statistical models**

All major steps are modularized in `vcf_process_func.py`,  
and orchestrated by the main driver script `filter_oriv.py`.
