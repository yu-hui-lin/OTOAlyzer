#!/usr/bin/env python3
#
# OTOAlyzer: CNV/SV and variant caller for OTOA gene and pseudogene from WGS
# Author: Yu-Hui Lin <yhlin.md05@nycu.edu.tw>
# OTOAlyzer is based on CNVPanelizer and Cyrius
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""
OTOAlyzer: OTOA gene CNV/SV and variant caller for WGS data

This tool determines:
1. Copy number of OTOA true gene and pseudogene (via cnvkit flat reference)
2. Structural variants (deletions, duplications, gene conversions)
3. De novo SNP/indel variants in OTOA true gene region

Pipeline:
=========
STEP 1: cnvkit (flat reference) → total_cn, true_gene_cn (OTOA_unique)
STEP 2: SNP counting at differentiating sites → per-site CN validation
STEP 3: Regional consensus analysis → SV detection
STEP 4: De novo variant calling → SNPs/indels across true gene
STEP 5: Quality filtering and output generation

CNV Calling Methods:
====================
- cnvkit (default): Uses flat reference with genome-wide bins for normalization.
                    No reference panel needed. Requires reference FASTA.
- CNVPanelizer (legacy): Uses panel of normal samples for normalization.
                         Requires reference BAM directory with 20+ diploid samples.

Output Files:
=============
- {prefix}.tsv  - Summary table (CN, SV, variant count)
- {prefix}.json - Detailed results with all variants
- {prefix}.vcf  - VCF format for variant integration
"""

import os
import sys
import argparse
import json
import logging
import datetime
from collections import namedtuple, OrderedDict
from typing import Dict, List, Tuple, Optional, Any
import glob

# Third-party imports
try:
    import pysam
except ImportError:
    print("ERROR: pysam is required. Install with: pip install pysam")
    sys.exit(1)

# =============================================================================
# Import OTOAlyzer modules
# =============================================================================

# SNP counting and CN calling (from Cyrius, adapted for OTOA)
from depth_calling.snp_count import (
    get_snp_position,
    get_supporting_reads,
    get_fraction,
)
from depth_calling.copy_number_call import (
    call_reg1_cn,
    process_raw_call_gc,
    process_raw_call_denovo,
)

# CNV calling - cnvkit (default) or CNVPanelizer (legacy)
from depth_calling.cnvkit_cn import get_cn_from_cnvkit

# Legacy CNVPanelizer integration (optional)
try:
    from depth_calling.panel_cn import get_cn_from_cnvpanelizer
    CNVPANELIZER_AVAILABLE = True
except ImportError:
    CNVPANELIZER_AVAILABLE = False

# OTOA-specific SV detection
from caller.otoa_sv import (
    get_otoa_sv_call,
    OtoaSvResult,
    get_sv_description,
    is_pathogenic_sv,
)

# De novo variant calling
from caller.otoa_variants import (
    call_variants_denovo,
    OtoaVariantCall,
    load_homology_sites,
    summarize_variants,
    variants_to_dict,
    filter_variants_by_quality,
    update_variants_with_phasing,
)

# Haplotype phasing
from depth_calling.haplotype import (
    phase_variants,
    summarize_phasing,
    detect_compound_het,
)


# =============================================================================
# Constants and Configuration
# =============================================================================

# Minimum read support thresholds
MIN_READ_SUPPORT = 5
MIN_CN_SITES_FRACTION = 0.5  # Minimum fraction of sites with valid CN calls

# Output data structures
OtoaCall = namedtuple(
    "OtoaCall",
    [
        "sample_id",
        "true_gene_cn",          # Copy number of OTOA true gene (from OTOA_unique)
        "pseudogene_cn",         # Copy number of OTOAP1 pseudogene
        "total_cn",              # Total CN (true + pseudo)
        "sv_call",               # Structural variant call
        "sv_confidence",         # Confidence score for SV
        "breakpoint",            # Breakpoint region if applicable
        "variants",              # List of called variants
        "phasing_summary",       # Haplotype phasing summary
        "compound_het",          # List of compound het variant pairs
        "filter_status",         # PASS or filter reason
        "raw_data",              # Raw supporting data for debugging
    ]
)

OtoaRawData = namedtuple(
    "OtoaRawData",
    [
        "cnv_method",             # 'cnvkit' or 'cnvpanelizer'
        "cnv_total_cn",           # Total CN from CNV caller
        "cnv_otoa_cn",            # OTOA CN from CNV caller
        "snp_true_gene_counts",
        "snp_pseudo_counts",
        "cn_per_site",
        "regional_consensus",
    ]
)


# =============================================================================
# Utility Functions
# =============================================================================

def open_alignment_file(bam_path: str, reference_fasta: Optional[str] = None) -> pysam.AlignmentFile:
    """
    Open a BAM/CRAM file with appropriate mode.
    
    Args:
        bam_path: Path to BAM or CRAM file
        reference_fasta: Path to reference FASTA (required for CRAM without embedded reference)
        
    Returns:
        pysam.AlignmentFile object
        
    Note:
        For CRAM files, the reference can be specified via:
        1. reference_fasta parameter
        2. REF_PATH environment variable
        3. Embedded in CRAM file header
    """
    if bam_path.endswith('.cram'):
        if reference_fasta:
            return pysam.AlignmentFile(bam_path, 'rc', reference_filename=reference_fasta)
        else:
            return pysam.AlignmentFile(bam_path, 'rc')
    else:
        return pysam.AlignmentFile(bam_path, 'rb')


def parse_bed_regions(bed_file: str) -> Dict[str, Dict]:
    """
    Parse BED file to extract region coordinates.
    
    Expected format:
    chr16  start  end  region_name
    
    Args:
        bed_file: Path to BED file
        
    Returns:
        Dictionary mapping region names to coordinates
    """
    regions = {}
    
    with open(bed_file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 4:
                chrom = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                name = fields[3]
                
                regions[name] = {
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'name': name,
                }
    
    return regions


def find_true_gene_region(region_dic: Dict) -> Optional[Dict]:
    """
    Find the true gene region from the BED file regions.
    
    Looks for region names containing 'OTOA_unique', 'OTOA_full', 'true', or 'OTOA'.
    
    Args:
        region_dic: Dictionary of regions from BED file
        
    Returns:
        Region dict or None
    """
    # Priority order for finding true gene region
    priority_patterns = [
        'OTOA_unique',      # Best - unique to true gene
        'OTOA_full',        # Full OTOA region
        'true_gene',        # Explicit true gene
        'true',             # Contains 'true'
        'OTOA',             # Just OTOA
    ]
    
    for pattern in priority_patterns:
        for name, region in region_dic.items():
            if pattern.lower() in name.lower():
                return region
    
    return None


def find_pseudogene_region(region_dic: Dict) -> Optional[Dict]:
    """
    Find the pseudogene region from the BED file regions.
    
    Args:
        region_dic: Dictionary of regions from BED file
        
    Returns:
        Region dict or None
    """
    patterns = ['OTOAP1_full', 'pseudogene', 'pseudo', 'OTOAP']
    
    for pattern in patterns:
        for name, region in region_dic.items():
            if pattern.lower() in name.lower():
                return region
    
    return None


# =============================================================================
# Resource Loading
# =============================================================================

def prepare_resources(datadir: str, genome: str, bed_file_path: str) -> Dict:
    """
    Load OTOA-specific resource files.
    
    Required files:
    - BED file defining genomic regions for CNVPanelizer
    - OTOA_SNP_{genome}.txt - differentiating SNPs between true gene and pseudogene
    
    Args:
        datadir: Path to data directory
        genome: Reference genome version (eg: '38' for GRCh38)
        bed_file_path: Path to BED file
        
    Returns:
        Dictionary containing parsed resource data
    """
    # SNP file for CN calling (REQUIRED)
    snp_file = os.path.join(datadir, f"OTOA_SNP_{genome}.txt")
    
    # Validate required files exist
    required_files = [bed_file_path, snp_file]
    for required_file in required_files:
        if not os.path.exists(required_file):
            raise FileNotFoundError(f"Required file not found: {required_file}")
    
    logging.info(f"Loading BED file: {bed_file_path}")
    logging.info(f"Loading SNP file: {snp_file}")
    
    # Parse region definitions from BED file
    region_dic = parse_bed_regions(bed_file_path)
    logging.info(f"Loaded {len(region_dic)} regions from BED file")
    
    # Parse SNP positions for CN calling
    snp_db = get_snp_position(snp_file, genome)
    logging.info(f"Loaded {len(snp_db.dsnp1)} differentiating SNP sites")
    
    # Load homology sites for variant filtering
    homology_sites = load_homology_sites(snp_file)
    logging.info(f"Loaded {len(homology_sites)} homology sites for variant filtering")
    
    # Find true gene and pseudogene regions
    true_gene_region = find_true_gene_region(region_dic)
    pseudogene_region = find_pseudogene_region(region_dic)
    
    if true_gene_region:
        logging.info(f"True gene region: {true_gene_region['name']} "
                     f"({true_gene_region['chrom']}:{true_gene_region['start']}-{true_gene_region['end']})")
    else:
        logging.warning("Could not identify true gene region in BED file")
    
    if pseudogene_region:
        logging.info(f"Pseudogene region: {pseudogene_region['name']} "
                     f"({pseudogene_region['chrom']}:{pseudogene_region['start']}-{pseudogene_region['end']})")
    
    return {
        "genome": genome,
        "region_dic": region_dic,
        "snp_db": snp_db,
        "homology_sites": homology_sites,
        "bed_file": bed_file_path,
        "true_gene_region": true_gene_region,
        "pseudogene_region": pseudogene_region,
    }


# =============================================================================
# Main Calling Logic
# =============================================================================

def call_sample(
    bam_file: str,
    resources: Dict,
    output_dir: str,
    reference_fasta: str,
    cnv_method: str = "cnvkit",
    reference_dir: Optional[str] = None,
    r_script_path: Optional[str] = None,
) -> OtoaCall:
    """
    Main function to call CNV/SV and variants for a single sample.

    Pipeline Steps:
    1. CNV calling (cnvkit or CNVPanelizer) → total_cn, true_gene_cn (OTOA_unique_cn)
    2. SNP counting at differentiating sites
    3. Per-site CN calling with Poisson model
    4. SV detection via regional consensus
    5. De novo variant calling
    6. Quality filtering

    Args:
        bam_file: Path to input BAM file
        resources: Dictionary of loaded resources
        output_dir: Output directory for intermediate files
        reference_fasta: Path to reference FASTA (REQUIRED for cnvkit, optional for CNVPanelizer)
        cnv_method: CNV calling method ('cnvkit' or 'cnvpanelizer')
        reference_dir: Path to reference BAM directory (only for CNVPanelizer)
        r_script_path: Path to CNVPanelizer R script (only for CNVPanelizer)

    Returns:
        OtoaCall named tuple with all results
    """
    sample_id = os.path.splitext(os.path.basename(bam_file))[0]
    logging.info(f"{'='*60}")
    logging.info(f"Processing sample: {sample_id}")
    logging.info(f"{'='*60}")
    
    # Open BAM/CRAM file
    try:
        bamfile = open_alignment_file(bam_file, reference_fasta=reference_fasta)
    except Exception as e:
        logging.error(f"Failed to open alignment file: {e}")
        if bam_file.endswith('.cram'):
            logging.error("For CRAM files, ensure reference FASTA is available via --reference or REF_PATH environment variable")
        return _create_failed_call(sample_id, f"Alignment_open_failed: {e}")
    
    # =========================================================================
    # STEP 1: Determine Copy Number via cnvkit or CNVPanelizer
    # =========================================================================
    logging.info(f"STEP 1: Determining copy number with {cnv_method}...")

    if cnv_method == "cnvkit":
        # cnvkit flat reference approach (default)
        if not reference_fasta:
            logging.error("Reference FASTA is REQUIRED for cnvkit")
            return _create_failed_call(sample_id, "Missing_reference_fasta")

        total_cn, otoa_cn = get_cn_from_cnvkit(
            bam_file=bam_file,
            bed_file_path=resources["bed_file"],
            reference_fasta=reference_fasta,
            output_dir=output_dir,
        )
        if total_cn is None or otoa_cn is None:
            logging.error(f"cnvkit failed for {sample_id}")
            return _create_failed_call(sample_id, "cnvkit_failed")

    elif cnv_method == "cnvpanelizer":
        # Legacy CNVPanelizer approach
        if not CNVPANELIZER_AVAILABLE:
            logging.error("CNVPanelizer module not available")
            return _create_failed_call(sample_id, "CNVPanelizer_not_available")

        if not reference_dir or not r_script_path:
            logging.error("reference_dir and r_script_path are required for CNVPanelizer")
            return _create_failed_call(sample_id, "Missing_CNVPanelizer_config")

        total_cn, otoa_cn = get_cn_from_cnvpanelizer(
            bam_file=bam_file,
            r_script_path=r_script_path,
            output_dir=output_dir,
            bed_file_path=resources["bed_file"],
            reference_dir_path=reference_dir,
            reference_fasta=reference_fasta,
        )

        if total_cn is None or otoa_cn is None:
            logging.error(f"CNVPanelizer failed for {sample_id}")
            return _create_failed_call(sample_id, "CNVPanelizer_failed")

    else:
        logging.error(f"Unknown CNV method: {cnv_method}")
        return _create_failed_call(sample_id, f"Unknown_cnv_method_{cnv_method}")
    
    # Calculate copy numbers
    # otoa_cn = OTOA_unique (true gene only)
    # total_cn = OTOA_full + OTOAP1_full (true + pseudo)
    true_gene_cn = otoa_cn
    pseudogene_cn = total_cn - otoa_cn
    
    logging.info(f"CNV calling results ({cnv_method}):")
    logging.info(f"  Total CN (OTOA+OTOAP1): {total_cn}")
    logging.info(f"  True gene CN (OTOA_unique): {true_gene_cn}")
    logging.info(f"  Pseudogene CN (calculated): {pseudogene_cn}")
    
    # Sanity check
    if pseudogene_cn < 0:
        logging.warning(f"Negative pseudogene CN ({pseudogene_cn}), setting to 0")
        pseudogene_cn = 0
    
    # =========================================================================
    # STEP 2: Count reads at differentiating SNP sites
    # =========================================================================
    
    # Open BAM/CRAM file
    try:
        bamfile = open_alignment_file(bam_file, reference_fasta=reference_fasta)
    except Exception as e:
        logging.error(f"Failed to open alignment file: {e}")
        if bam_file.endswith('.cram'):
            logging.error("For CRAM files, ensure reference FASTA is available via --reference or REF_PATH environment variable")
        return _create_failed_call(sample_id, f"Alignment_open_failed: {e}")
    logging.info("STEP 2: Counting reads at differentiating SNP sites...")
    
    snp_db = resources["snp_db"]
    
    # Get supporting reads for true gene (dsnp1) and pseudogene (dsnp2)
    snp_true_counts, snp_pseudo_counts = get_supporting_reads(
        bamfile,
        snp_db.dsnp1,  # True gene alleles
        snp_db.dsnp2,  # Pseudogene alleles
        snp_db.nchr,
        snp_db.dindex,
    )
    
    logging.info(f"Counted reads at {len(snp_true_counts)} SNP sites")
    
    # =========================================================================
    # STEP 3: Per-site CN calling with Poisson model
    # =========================================================================
    logging.info("STEP 3: Calling CN at each SNP site...")
    
    cn_call_per_site_raw = []
    for i in range(len(snp_true_counts)):
        site_cn = call_reg1_cn(
            full_cn=total_cn,
            count_reg1=snp_true_counts[i],
            count_reg2=snp_pseudo_counts[i],
            min_read=MIN_READ_SUPPORT,
        )
        cn_call_per_site_raw.append(site_cn)
    
    # Process raw calls to get consensus values
    cn_call_per_site = process_raw_call_gc(
        cn_prob=cn_call_per_site_raw,
        post_cutoff=0.8,
        keep_none=True,
    )
    
    # Count valid CN calls
    valid_cn_calls = [c for c in cn_call_per_site if c is not None]
    logging.info(f"Valid CN calls: {len(valid_cn_calls)}/{len(cn_call_per_site)} sites "
                 f"({100*len(valid_cn_calls)/len(cn_call_per_site):.1f}%)")
    
    # =========================================================================
    # STEP 4: Detect Structural Variants
    # =========================================================================
    logging.info("STEP 4: Detecting structural variants...")
    
    sv_result = get_otoa_sv_call(
        total_cn=total_cn,
        cn_call_per_site=cn_call_per_site,
        true_gene_cn=true_gene_cn,
        pseudogene_cn=pseudogene_cn,
    )
    
    logging.info(f"SV call: {sv_result.sv_type} (confidence: {sv_result.confidence:.2f})")
    if sv_result.breakpoint_region:
        logging.info(f"Breakpoint: {sv_result.breakpoint_region}")
    
    # =========================================================================
    # STEP 5: De Novo Variant Calling
    # =========================================================================
    logging.info("STEP 5: Calling variants de novo across true gene region...")
    
    true_gene_region = resources.get("true_gene_region")
    pseudogene_region = resources.get("pseudogene_region")
    
    # Initialize phasing variables
    phasing_summary = None
    compound_het_pairs = []

    if true_gene_region is None:
        logging.warning("Could not find true gene region - skipping variant calling")
        variants = []
    else:
        # Add true_gene_start to pseudogene region for contamination checking
        if pseudogene_region:
            pseudogene_region_copy = pseudogene_region.copy()
            pseudogene_region_copy['true_gene_start'] = true_gene_region['start']
        else:
            pseudogene_region_copy = None
        
        variants = call_variants_denovo(
            bamfile=bamfile,
            true_gene_region=true_gene_region,
            pseudogene_region=pseudogene_region_copy,
            homology_sites=resources.get("homology_sites"),
            true_gene_cn=true_gene_cn,
            reference_fasta=reference_fasta,
        )

        # Filter variants by quality
        pass_variants = [v for v in variants if v.filter_status == "PASS"]
        logging.info(f"Variants found: {len(variants)} total, {len(pass_variants)} PASS")

        # =====================================================================
        # STEP 5b: Haplotype Phasing (for multiple variants)
        # =====================================================================
        if len(pass_variants) >= 2:
            logging.info("STEP 5b: Performing haplotype phasing...")

            phasing_result = phase_variants(
                bamfile_handle=bamfile,
                variants=pass_variants,
                chromosome=true_gene_region['chrom'],
                snp_db=snp_db,
            )

            # Update variants with phasing information
            variants = update_variants_with_phasing(variants, phasing_result)

            # Get phasing summary
            phasing_summary = summarize_phasing(phasing_result, pass_variants)
            compound_het_pairs = phasing_summary.get("compound_het_variants", [])

            logging.info(f"Phasing: {phasing_summary['phase_sets']} phase sets, "
                        f"{phasing_summary['compound_het_pairs']} compound het pairs, "
                        f"{phasing_summary.get('pseudogene_reads_filtered', 0)} pseudogene reads filtered")
        else:
            logging.info("STEP 5b: Skipping phasing (fewer than 2 PASS variants)")
    
    # =========================================================================
    # STEP 6: Determine Filter Status
    # =========================================================================
    filter_status = determine_filter_status(
        true_gene_cn=true_gene_cn,
        pseudogene_cn=pseudogene_cn,
        sv_result=sv_result,
        cn_per_site=cn_call_per_site,
    )
    
    logging.info(f"Filter status: {filter_status}")
    
    # Compile raw data for debugging/output
    raw_data = OtoaRawData(
        cnv_method=cnv_method,
        cnv_total_cn=total_cn,
        cnv_otoa_cn=otoa_cn,
        snp_true_gene_counts=snp_true_counts,
        snp_pseudo_counts=snp_pseudo_counts,
        cn_per_site=cn_call_per_site,
        regional_consensus=sv_result.cn_pattern,
    )
    
    bamfile.close()
    
    return OtoaCall(
        sample_id=sample_id,
        true_gene_cn=true_gene_cn,
        pseudogene_cn=pseudogene_cn,
        total_cn=total_cn,
        sv_call=sv_result.sv_type,
        sv_confidence=sv_result.confidence,
        breakpoint=sv_result.breakpoint_region,
        variants=variants,
        phasing_summary=phasing_summary,
        compound_het=compound_het_pairs,
        filter_status=filter_status,
        raw_data=raw_data,
    )


def _create_failed_call(sample_id: str, reason: str) -> OtoaCall:
    """Create an OtoaCall for a failed sample."""
    return OtoaCall(
        sample_id=sample_id,
        true_gene_cn=None,
        pseudogene_cn=None,
        total_cn=None,
        sv_call=None,
        sv_confidence=None,
        breakpoint=None,
        variants=[],
        phasing_summary=None,
        compound_het=[],
        filter_status=reason,
        raw_data=None,
    )


def determine_filter_status(
    true_gene_cn: int,
    pseudogene_cn: int,
    sv_result: OtoaSvResult,
    cn_per_site: List[Optional[int]],
) -> str:
    """
    Determine QC filter status for the call.
    
    Args:
        true_gene_cn: True gene copy number
        pseudogene_cn: Pseudogene copy number
        sv_result: SV detection result
        cn_per_site: Per-site CN calls
        
    Returns:
        Filter status string (PASS or reason for filtering)
    """
    filters = []
    
    # Check for sufficient valid CN calls
    valid_cn_calls = [c for c in cn_per_site if c is not None]
    if len(valid_cn_calls) < len(cn_per_site) * MIN_CN_SITES_FRACTION:
        filters.append("LowConfidence_insufficient_sites")
    
    # Check SV confidence
    if sv_result.confidence is not None and sv_result.confidence < 0.6:
        filters.append("LowConfidence_SV")
    
    # Check for abnormally high CN
    if true_gene_cn is not None and pseudogene_cn is not None:
        if true_gene_cn + pseudogene_cn > 8:
            filters.append("LowQ_high_CN")
    
    # Check for negative pseudogene CN (shouldn't happen but sanity check)
    if pseudogene_cn is not None and pseudogene_cn < 0:
        filters.append("Invalid_negative_CN")
    
    return ";".join(filters) if filters else "PASS"


# =============================================================================
# Output Functions
# =============================================================================

def write_tsv_output(results: List[OtoaCall], output_path: str):
    """Write results to TSV file."""
    with open(output_path, 'w') as f:
        # Header
        f.write("Sample\tTrue_Gene_CN\tPseudogene_CN\tTotal_CN\t"
                "SV_Call\tSV_Confidence\tBreakpoint\t"
                "Num_Variants\tPASS_Variants\tVariants\tFilter\n")
        
        for call in results:
            # Count variants
            num_variants = len(call.variants) if call.variants else 0
            pass_variants = len([v for v in call.variants if v.filter_status == "PASS"]) if call.variants else 0
            conf_val = f"{call.sv_confidence:.2f}" if call.sv_confidence else "NA"
            # Format variant IDs
            variants_str = ";".join([v.variant_id for v in call.variants[:20]]) if call.variants else "None"
            if num_variants > 20:
                variants_str += f";...({num_variants-20} more)"
            
            f.write(f"{call.sample_id}\t"
                    f"{call.true_gene_cn if call.true_gene_cn is not None else 'NA'}\t"
                    f"{call.pseudogene_cn if call.pseudogene_cn is not None else 'NA'}\t"
                    f"{call.total_cn if call.total_cn is not None else 'NA'}\t"
                    f"{call.sv_call or 'NA'}\t"
                    f"{conf_val}\t"
                    f"{call.breakpoint or 'NA'}\t"
                    f"{num_variants}\t"
                    f"{pass_variants}\t"
                    f"{variants_str}\t"
                    f"{call.filter_status}\n")


def write_json_output(results: List[OtoaCall], output_path: str):
    """Write detailed results to JSON file."""
    output_data = OrderedDict()
    
    for call in results:
        sample_data = {
            "True_Gene_CN": call.true_gene_cn,
            "Pseudogene_CN": call.pseudogene_cn,
            "Total_CN": call.total_cn,
            "SV_Call": call.sv_call,
            "SV_Confidence": call.sv_confidence,
            "SV_Description": get_sv_description(call.sv_call) if call.sv_call else None,
            "Is_Pathogenic_SV": is_pathogenic_sv(call.sv_call) if call.sv_call else False,
            "Breakpoint": call.breakpoint,
            "Filter": call.filter_status,
            "Variants": [],
            "Variant_Summary": {},
        }
        
        # Add variant details
        if call.variants:
            sample_data["Variant_Summary"] = summarize_variants(call.variants)

            for var in call.variants:
                sample_data["Variants"].append({
                    "ID": var.variant_id,
                    "Chromosome": var.chromosome,
                    "Position": var.position,
                    "Ref": var.ref_base,
                    "Alt": var.alt_base,
                    "Type": var.variant_type,
                    "Allele_Fraction": var.allele_fraction,
                    "Zygosity": var.zygosity,
                    "Total_Depth": var.total_depth,
                    "Ref_Reads": var.ref_reads,
                    "Alt_Reads": var.alt_reads,
                    "Quality": var.quality_score,
                    "Filter": var.filter_status,
                    "In_Homology_Region": var.in_homology_region,
                    "Pseudogene_Evidence": var.pseudogene_evidence,
                    "Haplotype": var.haplotype,
                    "Phase_Set": var.phase_set,
                })

        # Add phasing summary
        if call.phasing_summary:
            sample_data["Phasing"] = {
                "Total_Pairs_Analyzed": call.phasing_summary.get("total_pairs", 0),
                "Cis_Pairs": call.phasing_summary.get("cis_pairs", 0),
                "Trans_Pairs": call.phasing_summary.get("trans_pairs", 0),
                "Unknown_Pairs": call.phasing_summary.get("unknown_pairs", 0),
                "Phase_Sets": call.phasing_summary.get("phase_sets", 0),
                "Compound_Het_Count": call.phasing_summary.get("compound_het_pairs", 0),
                "Pseudogene_Reads_Filtered": call.phasing_summary.get("pseudogene_reads_filtered", 0),
            }

        # Add compound heterozygous pairs
        if call.compound_het:
            sample_data["Compound_Heterozygous"] = [
                {"Variant1": pair[0], "Variant2": pair[1]}
                for pair in call.compound_het
            ]
        
        # Add raw data if available
        if call.raw_data:
            sample_data["Raw_Data"] = {
                "CNV_Method": call.raw_data.cnv_method,
                "CNV_Total_CN": call.raw_data.cnv_total_cn,
                "CNV_OTOA_CN": call.raw_data.cnv_otoa_cn,
                "Regional_Consensus": call.raw_data.regional_consensus,
            }
        
        output_data[call.sample_id] = sample_data
    
    with open(output_path, 'w') as f:
        json.dump(output_data, f, indent=2)


def write_vcf_output(results: List[OtoaCall], output_path: str, reference: str):
    """Write variants to VCF format."""
    with open(output_path, 'w') as f:
        # VCF Header
        f.write("##fileformat=VCFv4.2\n")
        f.write(f"##fileDate={datetime.datetime.now().strftime('%Y%m%d')}\n")
        f.write("##source=OTOAlyzer\n")
        f.write(f"##reference={reference}\n")
        f.write('##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant type (SNV, INS, DEL)">\n')
        f.write('##INFO=<ID=HOMOLOGY,Number=0,Type=Flag,Description="Variant in homology region">\n')
        f.write('##INFO=<ID=PSEUDO,Number=0,Type=Flag,Description="Possible pseudogene contamination">\n')
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write('##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele fraction">\n')
        f.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n')
        f.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths (ref,alt)">\n')
        f.write('##FILTER=<ID=LowDepth,Description="Total depth below threshold">\n')
        f.write('##FILTER=<ID=LowAltReads,Description="Alt read count below threshold">\n')
        f.write('##FILTER=<ID=StrandBias,Description="Significant strand bias">\n')
        f.write('##FILTER=<ID=HomologyRegion,Description="Position in high-homology region">\n')
        f.write('##FILTER=<ID=PseudogeneContamination,Description="Possible pseudogene contamination">\n')
        f.write('##FILTER=<ID=BorderlineAF,Description="Borderline allele fraction">\n')
        
        # Sample columns
        sample_ids = [r.sample_id for r in results]
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + 
                "\t".join(sample_ids) + "\n")
        
        # Collect all unique variants across samples
        all_variants = {}
        for call in results:
            if call.variants:
                for var in call.variants:
                    key = (var.chromosome, var.position, var.ref_base, var.alt_base)
                    if key not in all_variants:
                        all_variants[key] = var
        
        # Write variant lines (sorted by position)
        for (chrom, pos, ref, alt), var in sorted(all_variants.items(), key=lambda x: (x[0][0], x[0][1])):
            # Build INFO field
            info_parts = [f"TYPE={var.variant_type}"]
            if var.in_homology_region:
                info_parts.append("HOMOLOGY")
            if var.pseudogene_evidence:
                info_parts.append("PSEUDO")
            info_str = ";".join(info_parts) if info_parts else "."
            
            # Get genotypes for each sample
            sample_genotypes = []
            for call in results:
                sample_var = None
                if call.variants:
                    for v in call.variants:
                        if (v.chromosome, v.position, v.ref_base, v.alt_base) == (chrom, pos, ref, alt):
                            sample_var = v
                            break
                
                if sample_var:
                    # Determine GT from zygosity
                    if sample_var.zygosity == "HOM":
                        gt = "1/1"
                    elif "HET" in sample_var.zygosity:
                        gt = "0/1"
                    elif sample_var.zygosity == "HEMI":
                        gt = "1"
                    else:
                        gt = "./."
                    
                    ad = f"{sample_var.ref_reads},{sample_var.alt_reads}"
                    sample_genotypes.append(f"{gt}:{sample_var.allele_fraction:.3f}:{sample_var.total_depth}:{ad}")
                else:
                    sample_genotypes.append("0/0:0.000:0:0,0")
            
            # Write line
            qual = f"{var.quality_score:.1f}"
            filt = var.filter_status if var.filter_status != "PASS" else "PASS"
            
            f.write(f"{chrom}\t{pos}\t{var.variant_id}\t{ref}\t{alt}\t{qual}\t"
                    f"{filt}\t{info_str}\tGT:AF:DP:AD\t" + "\t".join(sample_genotypes) + "\n")


# =============================================================================
# Main Entry Point
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="OTOAlyzer: CNV/SV and de novo variant caller for OTOA gene from WGS",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Using cnvkit (default, recommended)
  python otoa_caller.py --manifest samples.txt --genome 38 --prefix output \\
                        --outDir results/ --reference hg38.fa

  # Using CNVPanelizer (legacy, requires reference panel)
  python otoa_caller.py --manifest samples.txt --genome 38 --prefix output \\
                        --outDir results/ --cnv-method cnvpanelizer --refDir ref_bams/

  # With custom data directory
  python otoa_caller.py --manifest samples.txt --genome 38 --prefix output \\
                        --outDir results/ --reference hg38.fa --dataDir /path/to/data/
        """
    )
    parser.add_argument(
        "--manifest",
        required=True,
        help="Path to manifest file listing BAM/CRAM files (one per line)",
    )
    parser.add_argument(
        "--genome",
        required=True,
        choices=["38","37"],
        help="Reference genome version (currently only GRCh37 and GRCh38 supported)",
    )
    parser.add_argument(
        "--prefix",
        required=True,
        help="Output file prefix",
    )
    parser.add_argument(
        "--outDir",
        required=True,
        help="Output directory",
    )
    parser.add_argument(
        "--cnv_method",
        choices=["cnvkit", "cnvpanelizer"],
        default="cnvkit",
        help="CNV calling method. 'cnvkit' (default) uses flat reference with genome-wide "
             "normalization (requires --reference). 'cnvpanelizer' uses panel of normal samples "
             "(requires --refDir with 20+ diploid BAMs).",
    )
    parser.add_argument(
        "--refDir",
        default=None,
        help="Directory containing reference BAM files for CNVPanelizer (20+ diploid samples). "
             "Only required when --cnv-method=cnvpanelizer. Default: ./ref directory.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of threads (default: 1) [not yet implemented]",
    )
    parser.add_argument(
        "--bed",
        help="Path to BED file (auto-detected from data directory if not specified)",
    )
    parser.add_argument(
        "--dataDir",
        default=None,
        help="Path to data directory containing OTOA_SNP_genome.txt (default: ./data)",
    )
    parser.add_argument(
        "--reference",
        help="Path to reference FASTA. REQUIRED when using cnvkit (default) for GC correction, "
             "and REQUIRED for CRAM input files. Recommended for accurate variant calling.",
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose logging",
    )
    
    args = parser.parse_args()
    
    # Setup logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    
    logging.info("="*60)
    logging.info("OTOAlyzer: OTOA CNV/SV and Variant Caller")
    logging.info("="*60)

    # Get CNV method
    cnv_method = getattr(args, 'cnv_method', 'cnvkit')
    logging.info(f"CNV calling method: {cnv_method}")

    # Validate method-specific requirements
    if cnv_method == "cnvkit":
        if not args.reference:
            logging.error("Reference FASTA (--reference) is REQUIRED for cnvkit")
            logging.error("The reference is needed for GC content correction in cnvkit")
            sys.exit(1)
        if not os.path.exists(args.reference):
            logging.error(f"Reference FASTA not found: {args.reference}")
            sys.exit(1)
        logging.info(f"Reference FASTA: {args.reference}")

    # Determine paths
    script_dir = os.path.dirname(os.path.realpath(__file__))
    data_dir = args.dataDir or os.path.join(script_dir, "data")

    logging.info(f"Data directory: {data_dir}")

    # Reference directory (only needed for CNVPanelizer)
    ref_dir = None
    r_script_path = None

    if cnv_method == "cnvpanelizer":
        if args.refDir:
            ref_dir = args.refDir
        else:
            ref_dir = os.path.join(script_dir, "ref")
            logging.info(f"Using default reference directory: {ref_dir}")

        logging.info(f"Reference BAM directory: {ref_dir}")

        # R script path for CNVPanelizer
        r_script_path = os.path.join(script_dir, "depth_calling", "run_CNVPanelizer.R")
        if not os.path.exists(r_script_path):
            logging.error(f"CNVPanelizer R script not found: {r_script_path}")
            sys.exit(1)

        # Validate reference directory
        if not os.path.isdir(ref_dir):
            logging.error(f"Reference directory not found: {ref_dir}")
            logging.error("Please create the ref/ directory and add 20+ diploid reference BAMs, "
                          "or specify --refDir")
            sys.exit(1)

        ref_bams = glob.glob(os.path.join(ref_dir, "*.bam")) + glob.glob(os.path.join(ref_dir, "*.cram"))
        logging.info(f"Found {len(ref_bams)} reference alignment files (BAM/CRAM)")
        if len(ref_bams) < 20:
            logging.warning(f"CNVPanelizer works best with 20+ reference samples (found {len(ref_bams)})")
    
    # Auto-detect BED file if not specified
    if args.bed:
        bed_file_path = args.bed
    else:
        # Look for BED files in data directory
        bed_patterns = [
            os.path.join(data_dir, "*region*.bed"),
            os.path.join(data_dir, "*OTOA*.bed"),
            os.path.join(data_dir, "*.bed"),
        ]
        bed_files = []
        for pattern in bed_patterns:
            bed_files.extend(glob.glob(pattern))
        
        if not bed_files:
            logging.error(f"No BED file found in {data_dir}. Please specify --bed")
            sys.exit(1)
        
        bed_file_path = bed_files[0]
        logging.info(f"Auto-detected BED file: {bed_file_path}")

    # Load resources
    logging.info("Loading resources...")
    try:
        resources = prepare_resources(data_dir, args.genome, bed_file_path)
    except FileNotFoundError as e:
        logging.error(str(e))
        sys.exit(1)
    
    # Create output directory
    os.makedirs(args.outDir, exist_ok=True)
    
    # Read manifest
    if not os.path.exists(args.manifest):
        logging.error(f"Manifest file not found: {args.manifest}")
        sys.exit(1)
    
    with open(args.manifest, 'r') as f:
        bam_files = [line.strip() for line in f if line.strip() and not line.startswith('#')]

    # Validate CRAM files require reference FASTA
    cram_files = [f for f in bam_files if f.endswith('.cram')]
    if cram_files and not args.reference:
        logging.error("Reference FASTA (--reference) is REQUIRED for CRAM files")
        logging.error(f"Found {len(cram_files)} CRAM file(s) in manifest:")
        for cf in cram_files[:5]:  # Show first 5
            logging.error(f"  - {os.path.basename(cf)}")
        if len(cram_files) > 5:
            logging.error(f"  ... and {len(cram_files) - 5} more")
        sys.exit(1)

    logging.info(f"Processing {len(bam_files)} sample(s)...")
    
    # Process each sample
    results = []
    for i, bam_file in enumerate(bam_files, 1):
        logging.info(f"\n[{i}/{len(bam_files)}] Processing: {os.path.basename(bam_file)}")
        
        if not os.path.exists(bam_file):
            logging.error(f"BAM file not found: {bam_file}")
            continue
        
        try:
            call = call_sample(
                bam_file=bam_file,
                resources=resources,
                output_dir=args.outDir,
                reference_fasta=args.reference,
                cnv_method=cnv_method,
                reference_dir=ref_dir,
                r_script_path=r_script_path,
            )
            results.append(call)
            
            # Log summary for this sample
            logging.info(f"Sample {call.sample_id} completed:")
            logging.info(f"  CN: OTOA={call.true_gene_cn}, OTOAP1={call.pseudogene_cn}")
            conf_str = f"{call.sv_confidence:.2f}" if call.sv_confidence else "NA"
            logging.info(f" SV: {call.sv_call} (conf={conf_str})")
            logging.info(f"  Variants: {len(call.variants)} total")
            
        except Exception as e:
            logging.error(f"Error processing {bam_file}: {e}")
            import traceback
            traceback.print_exc()
    
    # Write outputs
    if results:
        tsv_path = os.path.join(args.outDir, f"{args.prefix}.tsv")
        json_path = os.path.join(args.outDir, f"{args.prefix}.json")
        vcf_path = os.path.join(args.outDir, f"{args.prefix}.vcf")
        
        logging.info("\nWriting output files...")
        write_tsv_output(results, tsv_path)
        write_json_output(results, json_path)
        write_vcf_output(results, vcf_path, str(args.genome))
        
        logging.info(f"\nResults written to:")
        logging.info(f"  TSV:  {tsv_path}")
        logging.info(f"  JSON: {json_path}")
        logging.info(f"  VCF:  {vcf_path}")
    else:
        logging.warning("No samples were successfully processed!")
    
    logging.info("\n" + "="*60)
    logging.info("OTOAlyzer completed!")
    logging.info("="*60)


if __name__ == "__main__":
    main()
