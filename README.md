# OTOAlyzer

An accurate CNV/SV and variant caller for OTOA gene and pseudogene from WGS data.

---

## Overview

OTOAlyzer determines:
1. **Copy Number** - True gene and pseudogene copy numbers
2. **Structural Variants** - Deletions, duplications, gene conversions, hybrids
3. **De Novo Variant Detection** - Scans entire OTOA true gene region for SNPs/indels

---

## Directory Structure

```
OTOAlyzer/
├── otoa_caller.py              # Main entry point
├── caller/
│   ├── __init__.py
│   ├── otoa_sv.py              # SV detection (regional consensus)
│   └── otoa_variants.py        # De novo variant calling
├── depth_calling/
│   ├── __init__.py
│   ├── panel_cn.py             # CNVPanelizer integration
│   ├── run_CNVPanelizer.R      # R script for CNVPanelizer
│   ├── copy_number_call.py     # Poisson CN calling
│   └── snp_count.py            # Read counting at SNP sites
├── data/
│   ├── OTOA_region_38.bed      # BED file for CNVPanelizer regions
│   └── OTOA_SNP_38.txt         # 160 differentiating SNPs for CN calling
└──  ref/                        # Reference BAMs for CNVPanelizer (20+ diploid samples)
    ├── sample01.bam
    ├── sample01.bam.bai
    └── ... (20+ samples)
```

---
