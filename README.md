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

## Setup Instructions


### Step 3: Prepare 20+ diploid samples stored in OTOAlyzer/ref/ to serve as copy number estimation
- Same sequencing platform as test samples
- Same alignment pipeline (GRCh38/hg38)
- All must have index files (.bai or .crai)
- No known OTOA CNVs or SVs

---

## Usage

OTOAlyzer supports both BAM and CRAM files. 
For CRAM files: --reference /path/to/GRCh38.fa is needed.

```bash
python3 otoa_caller.py \
    --manifest sample.manifest \   # Required: BAM/CRAM file list
    --genome 38 \                  # Required: Reference genome (38 only)
    --prefix sample_output \       # Required: Output file prefix
    --outDir results/ \            # Required: Output directory
    --refDir ref/ \                # Optional: Directory of reference diploid BAMs/CRAMs (default: ./ref)
    --dataDir data/ \              # Optional: Data files (default: ./data)
    --bed custom.bed \             # Optional: Custom BED file
    --reference hg38.fa \          # Optional: Reference FASTA
    --verbose                      # Optional: Verbose logging
```

### Input Files
**Manifest file:** Text file with one BAM/CRAM path per line
```
/path/to/sample1.bam
/path/to/sample2.cram
```

### Output Format

### TSV Output

```
Sample	True_Gene_CN	Pseudogene_CN	Total_CN	SV_Call	Num_Variants	Variants	Filter
NA12878	2	2	4	normal	3	chr16:21680123:A>G;chr16:21681234:C>T;chr16:21682000:G>A	PASS
NA12879	1	2	3	true_gene_del_het	0	None	PASS
```

### JSON Output

```json
{
    "NA12878": {
        "True_Gene_CN": 2,
        "Pseudogene_CN": 2,
        "Total_CN": 4,
        "SV_Call": "normal",
        "SV_Confidence": 0.95,
        "Filter": "PASS",
        "Variants": [
            {
                "ID": "chr16:21680123:A>G",
                "Position": 21680123,
                "Ref": "A",
                "Alt": "G",
                "Type": "SNV",
                "Allele_Fraction": 0.48,
                "Zygosity": "HET",
                "Total_Depth": 150,
                "Alt_Reads": 72,
                "Quality": 245.3,
                "Filter": "PASS",
                "In_Homology_Region": false,
                "Pseudogene_Evidence": null
            }
        ]
    }
}
```

### VCF Output

Standard VCF 4.2 format with:
- Per-sample genotypes (GT), allele fractions (AF), depths (DP), allelic depths (AD)
- INFO fields: TYPE (SNV/INS/DEL), HOMOLOGY flag, PSEUDO flag
- FILTER annotations for quality issues





