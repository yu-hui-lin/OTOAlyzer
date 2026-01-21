#!/usr/bin/env python3
#
# CNV and SV detection for OTOA and its pseudogene from WGS
# Author: Yu-Hui Lin <yhlin.md05@nycu.edu.tw>
# OTOAlyzer is based on CNVPanelizer and Cyrius
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""
Structural Variant Detection for OTOA

This module detects SVs in the OTOA gene region by analyzing
copy number patterns across the gene using regional consensus,
adapted from Cyrius's approach for CYP2D6.

OTOA Gene Structure (based on OTOA_SNP_38.txt annotations):
The homologous region between OTOA and OTOAP1 spans exon21-exon29.

SNP Distribution in OTOA_SNP_38.txt (157 SNPs total):
- upstreamparalog: 12 SNPs (indices 0-11)
- upstreamparalog_150bp: 2 SNPs (indices 12-13)
- exon21: 1 SNP (index 14)
- intron21: 63 SNPs (indices 15-77)
- exon22: 3 SNPs (indices 78-80)
- intron22: 17 SNPs (indices 81-97)
- intron23: 4 SNPs (indices 98-101)
- intron24: 10 SNPs (indices 102-111)
- intron26: 1 SNP (index 112)
- intron27: 12 SNPs (indices 113-124)
- intron28: 27 SNPs (indices 125-151)
- exon29: 3 SNPs (indices 152-154)
- down_exon29: 2 SNPs (indices 155-156)
Total: 157 SNPs

Regional Division (4 regions, like CYP2D6):
1. upstream_region:  upstreamparalog + upstreamparalog_150bp (14 SNPs, indices 0-13)
2. 5prime_region:    exon21 + intron21 (64 SNPs, indices 14-77)
3. middle_region:    exon22 + intron22 + intron23 + intron24 (34 SNPs, indices 78-111)
4. 3prime_region:    intron26 + intron27 + intron28 + exon29 + down_exon29 (45 SNPs, indices 112-156)

SV Types Detected:
- Deletions (heterozygous/homozygous) of true gene
- Deletions of pseudogene
- Duplications of true gene
- Gene conversions (partial CN changes)
- Hybrid genes
"""

import logging
from collections import namedtuple, Counter
from typing import List, Optional, Tuple, Dict

# =============================================================================
# Data Structures
# =============================================================================

OtoaSvResult = namedtuple(
    "OtoaSvResult",
    [
        "sv_type",           # Type of SV detected
        "confidence",        # Confidence score (0-1)
        "breakpoint_region", # Approximate breakpoint location
        "cn_pattern",        # CN pattern across regions
        "details",           # Additional details dict
    ]
)

OtoaRegionalCN = namedtuple(
    "OtoaRegionalCN",
    [
        "upstream",      # upstream paralog region (indices 0-13, 14 SNPs)
        "region_5prime", # exon21 + intron21 (indices 14-77, 64 SNPs)
        "region_middle", # exon22 through intron24 (indices 78-111, 34 SNPs)
        "region_3prime", # intron26 through exon29 (indices 112-156, 45 SNPs)
    ]
)

# =============================================================================
# Constants for OTOA Regional Analysis
# Based on OTOA_SNP_38.txt annotation column
# =============================================================================

# Region boundaries (SNP indices)
UPSTREAM_START = 0
UPSTREAM_END = 14           # indices 0-13: upstreamparalog + upstreamparalog_150bp (14 SNPs)

REGION_5PRIME_START = 14
REGION_5PRIME_END = 78      # indices 14-77: exon21 + intron21 (64 SNPs)

REGION_MIDDLE_START = 78
REGION_MIDDLE_END = 112     # indices 78-111: exon22 + intron22 + intron23 + intron24 (34 SNPs)

REGION_3PRIME_START = 112
REGION_3PRIME_END = 157     # indices 112-156: intron26 + intron27 + intron28 + exon29 + down_exon29 (45 SNPs)

# Minimum sites for reliable consensus (approximately 50% of region size)
UPSTREAM_SITES_MIN = 7          # ~50% of 14 SNPs
REGION_5PRIME_SITES_MIN = 32    # ~50% of 64 SNPs
REGION_MIDDLE_SITES_MIN = 17    # ~50% of 34 SNPs
REGION_3PRIME_SITES_MIN = 22    # ~50% of 45 SNPs

# Looser threshold when primary consensus fails
CONSENSUS_MIN_LOOSE = 10


# =============================================================================
# Main SV Detection
# =============================================================================

def get_otoa_sv_call(
    total_cn: int,
    cn_call_per_site: List[Optional[int]],
    true_gene_cn: int,
    pseudogene_cn: int,
) -> OtoaSvResult:
    """
    Determine SV type based on regional CN patterns.
    
    This is the main entry point for SV detection. It analyzes
    CN patterns across the OTOA gene regions to detect SVs.
    
    Args:
        total_cn: Total copy number (true gene + pseudogene)
        cn_call_per_site: Per-site CN calls from Poisson model
        true_gene_cn: CN of true OTOA gene (from CNVPanelizer OTOA_unique)
        pseudogene_cn: CN of pseudogene (total_cn - true_gene_cn)
        
    Returns:
        OtoaSvResult with SV classification
    """
    logging.info(f"SV detection: total_cn={total_cn}, true_gene={true_gene_cn}, pseudo={pseudogene_cn}")
    
    # Handle None values
    if true_gene_cn is None:
        return OtoaSvResult(
            sv_type="unknown",
            confidence=0.0,
            breakpoint_region=None,
            cn_pattern=None,
            details={"error": "Missing true gene CN"},
        )
    
    # Get regional consensus from per-site CN calls
    consensus = get_regional_consensus(cn_call_per_site, total_cn)
    
    logging.info(f"Regional consensus: upstream={consensus.upstream}, "
                 f"5prime={consensus.region_5prime}, middle={consensus.region_middle}, "
                 f"3prime={consensus.region_3prime}")
    
    # Classify SV based on CN pattern
    sv_type, confidence, change_points = classify_sv(
        total_cn=total_cn,
        true_gene_cn=true_gene_cn,
        pseudogene_cn=pseudogene_cn,
        consensus=consensus,
    )
    
    # Detect breakpoint if SV present
    breakpoint = None
    if sv_type not in [None, "normal", "cn2"]:
        breakpoint = detect_breakpoint(consensus, change_points)
    
    return OtoaSvResult(
        sv_type=sv_type,
        confidence=confidence,
        breakpoint_region=breakpoint,
        cn_pattern={
            "upstream": consensus.upstream,
            "5prime": consensus.region_5prime,
            "middle": consensus.region_middle,
            "3prime": consensus.region_3prime,
        },
        details={
            "total_cn": total_cn,
            "true_gene_cn": true_gene_cn,
            "pseudogene_cn": pseudogene_cn,
            "change_points": change_points,
        },
    )


def get_regional_consensus(
    cn_call_per_site: List[Optional[int]],
    total_cn: int,
) -> OtoaRegionalCN:
    """
    Calculate consensus CN for each region of the OTOA gene.
    
    Similar to Cyrius's approach for CYP2D6:
    - Extract CN calls for each region
    - Find most common CN value
    - Require minimum number of supporting sites
    
    Args:
        cn_call_per_site: Per-site CN calls from Poisson model
        total_cn: Total CN for reference
        
    Returns:
        OtoaRegionalCN named tuple with consensus CN for each region
    """
    # Handle empty or short input
    if not cn_call_per_site:
        return OtoaRegionalCN(None, None, None, None)
    
    # Extract sites for each region based on OTOA_SNP_38.txt structure
    upstream_sites = [
        a for a in cn_call_per_site[UPSTREAM_START:UPSTREAM_END]
        if a is not None
    ]
    
    region_5prime_sites = [
        a for a in cn_call_per_site[REGION_5PRIME_START:REGION_5PRIME_END]
        if a is not None
    ]
    
    region_middle_sites = [
        a for a in cn_call_per_site[REGION_MIDDLE_START:REGION_MIDDLE_END]
        if a is not None
    ]
    
    region_3prime_sites = [
        a for a in cn_call_per_site[REGION_3PRIME_START:REGION_3PRIME_END]
        if a is not None
    ]
    
    # Calculate consensus for each region
    upstream_consensus = _get_region_consensus(
        upstream_sites, UPSTREAM_SITES_MIN, total_cn
    )
    
    consensus_5prime = _get_region_consensus(
        region_5prime_sites, REGION_5PRIME_SITES_MIN, total_cn
    )
    
    consensus_middle = _get_region_consensus(
        region_middle_sites, REGION_MIDDLE_SITES_MIN, total_cn
    )
    
    consensus_3prime = _get_region_consensus(
        region_3prime_sites, REGION_3PRIME_SITES_MIN, total_cn
    )
    
    # Fill in missing consensus values from adjacent regions (like Cyrius)
    # This helps when one region has insufficient data
    if consensus_5prime is None and consensus_middle is not None:
        consensus_5prime = consensus_middle
    elif consensus_middle is None and consensus_5prime is not None:
        consensus_middle = consensus_5prime
    
    if consensus_middle is None and consensus_3prime is not None:
        consensus_middle = consensus_3prime
    elif consensus_3prime is None and consensus_middle is not None:
        consensus_3prime = consensus_middle
    
    return OtoaRegionalCN(
        upstream=upstream_consensus,
        region_5prime=consensus_5prime,
        region_middle=consensus_middle,
        region_3prime=consensus_3prime,
    )


def _get_region_consensus(
    site_calls: List[int],
    min_sites: int,
    total_cn: int,
) -> Optional[int]:
    """
    Get consensus CN from a list of site calls.
    
    Args:
        site_calls: List of CN values at each site
        min_sites: Minimum sites required for consensus
        total_cn: Total CN for fallback logic
        
    Returns:
        Consensus CN or None if insufficient data
    """
    if not site_calls:
        return None
    
    counter = sorted(
        Counter(site_calls).items(),
        key=lambda kv: kv[1],
        reverse=True
    )
    
    if counter:
        # Primary: require min_sites for consensus
        if counter[0][1] >= min_sites:
            return counter[0][0]
        
        # Fallback: looser threshold
        if counter[0][1] >= CONSENSUS_MIN_LOOSE:
            return counter[0][0]
        
        # Last resort: if expected CN (total_cn - 2 for pseudogene = 2) has reasonable support
        expected_otoa_cn = total_cn - 2
        if site_calls.count(expected_otoa_cn) >= CONSENSUS_MIN_LOOSE:
            return expected_otoa_cn
    
    return None


# =============================================================================
# SV Classification
# =============================================================================

def classify_sv(
    total_cn: int,
    true_gene_cn: int,
    pseudogene_cn: int,
    consensus: OtoaRegionalCN,
) -> Tuple[str, float, List[str]]:
    """
    Classify SV type based on CN patterns.
    
    Similar to Cyrius's get_cnvtag() function:
    - Analyze CN changes between regions
    - Map CN patterns to SV types
    
    Args:
        total_cn: Total CN
        true_gene_cn: True gene CN (from CNVPanelizer OTOA_unique)
        pseudogene_cn: Pseudogene CN
        consensus: Regional consensus values
        
    Returns:
        Tuple of (sv_type, confidence, change_points)
    """
    change_points = []
    
    # --- Check for normal diploid state ---
    # Normal: true gene = 2, pseudogene = 2, all regions = 2
    if true_gene_cn == 2 and pseudogene_cn == 2:
        all_regions_normal = all(
            c is None or c == 2
            for c in [consensus.upstream, consensus.region_5prime,
                      consensus.region_middle, consensus.region_3prime]
        )
        if all_regions_normal:
            return ("normal", 0.95, [])
    
    # --- Simple whole-gene events (based on CNVPanelizer results) ---
    
    # Heterozygous deletion of true gene (CN=1)
    if true_gene_cn == 1 and pseudogene_cn == 2:
        return ("true_gene_del_het", 0.90, ["whole_gene_deletion"])
    
    # Homozygous deletion of true gene (CN=0)
    if true_gene_cn == 0 and pseudogene_cn == 2:
        return ("true_gene_del_hom", 0.95, ["whole_gene_deletion_hom"])
    
    # Duplication of true gene (CN=3)
    if true_gene_cn == 3 and pseudogene_cn == 2:
        return ("true_gene_dup", 0.85, ["whole_gene_duplication"])
    
    # Triplication or higher of true gene
    if true_gene_cn >= 4 and pseudogene_cn == 2:
        return (f"true_gene_cn{true_gene_cn}", 0.80, ["whole_gene_multiplication"])
    
    # Heterozygous deletion of pseudogene
    if true_gene_cn == 2 and pseudogene_cn == 1:
        return ("pseudogene_del_het", 0.85, ["pseudogene_deletion"])
    
    # Homozygous deletion of pseudogene
    if true_gene_cn == 2 and pseudogene_cn == 0:
        return ("pseudogene_del_hom", 0.90, ["pseudogene_deletion_hom"])
    
    # --- Regional CN changes (potential gene conversion or partial events) ---
    # Analyze CN transitions between adjacent regions
    
    regions = [
        ("upstream", consensus.upstream),
        ("5prime", consensus.region_5prime),
        ("middle", consensus.region_middle),
        ("3prime", consensus.region_3prime),
    ]
    
    # Detect CN changes between regions
    for i in range(len(regions) - 1):
        region1_name, region1_cn = regions[i]
        region2_name, region2_cn = regions[i + 1]
        
        if region1_cn is not None and region2_cn is not None:
            if region1_cn > region2_cn:
                # CN decreases: potential partial deletion or gene conversion to pseudo
                change_points.append(f"del_{region1_name}_to_{region2_name}")
            elif region1_cn < region2_cn:
                # CN increases: potential partial duplication or gene conversion from pseudo
                change_points.append(f"dup_{region1_name}_to_{region2_name}")
    
    # --- Classify based on change points ---
    
    if change_points:
        # Gene conversion: CN changes within the gene
        if len(change_points) == 1:
            if "del" in change_points[0]:
                return ("partial_deletion", 0.75, change_points)
            else:
                return ("partial_duplication", 0.75, change_points)
        else:
            # Multiple change points: complex event
            return ("gene_conversion", 0.70, change_points)
    
    # --- Handle unusual patterns ---
    
    # Both true gene and pseudogene affected
    if true_gene_cn != 2 and pseudogene_cn != 2:
        if true_gene_cn < 2 and pseudogene_cn < 2:
            return ("combined_deletion", 0.70, ["both_deleted"])
        elif true_gene_cn > 2 and pseudogene_cn > 2:
            return ("combined_duplication", 0.65, ["both_duplicated"])
        else:
            return ("hybrid", 0.60, ["complex_rearrangement"])
    
    # Check for gene conversion when CNVPanelizer shows normal but regions vary
    if true_gene_cn == 2 and pseudogene_cn == 2:
        if has_regional_variation(consensus):
            return ("gene_conversion", 0.70, ["regional_variation"])
    
    # Cannot determine
    logging.warning(f"Could not classify SV: true_cn={true_gene_cn}, pseudo_cn={pseudogene_cn}, "
                    f"consensus={consensus}")
    return (None, 0.0, [])


def has_regional_variation(consensus: OtoaRegionalCN) -> bool:
    """
    Check if there's variation between gene regions suggesting gene conversion.
    """
    regions = [consensus.upstream, consensus.region_5prime, 
               consensus.region_middle, consensus.region_3prime]
    valid_regions = [r for r in regions if r is not None]
    
    if len(valid_regions) < 2:
        return False
    
    return len(set(valid_regions)) > 1


# =============================================================================
# Breakpoint Detection
# =============================================================================

def detect_breakpoint(
    consensus: OtoaRegionalCN,
    change_points: List[str],
) -> Optional[str]:
    """
    Detect approximate breakpoint location based on regional CN changes.
    
    Args:
        consensus: Regional consensus values
        change_points: List of detected change points
        
    Returns:
        String describing breakpoint region, or None
    """
    if not change_points:
        return None
    
    # Map change points to approximate genomic regions based on OTOA structure
    breakpoint_map = {
        "del_upstream_to_5prime": "upstream_exon21_boundary",
        "dup_upstream_to_5prime": "upstream_exon21_boundary",
        "del_5prime_to_middle": "intron21_exon22_boundary",
        "dup_5prime_to_middle": "intron21_exon22_boundary",
        "del_middle_to_3prime": "intron24_intron26_boundary",
        "dup_middle_to_3prime": "intron24_intron26_boundary",
    }
    
    breakpoints = []
    for cp in change_points:
        if cp in breakpoint_map:
            breakpoints.append(breakpoint_map[cp])
        else:
            breakpoints.append(cp)
    
    if breakpoints:
        return ";".join(breakpoints)
    
    return None


def get_detailed_breakpoint(
    cn_call_per_site: List[Optional[int]],
    snp_positions: List[int],
    region_start_idx: int,
    region_end_idx: int,
) -> Optional[Dict]:
    """
    Find exact SNP index where CN changes within a region.
    
    Args:
        cn_call_per_site: Per-site CN calls
        snp_positions: List of genomic positions for each SNP
        region_start_idx: Start index of region
        region_end_idx: End index of region
        
    Returns:
        Dict with breakpoint info, or None
    """
    region_calls = cn_call_per_site[region_start_idx:region_end_idx]
    
    prev_cn = None
    for i, cn in enumerate(region_calls):
        if cn is not None:
            if prev_cn is not None and cn != prev_cn:
                global_idx = region_start_idx + i
                if snp_positions and len(snp_positions) > global_idx:
                    return {
                        "snp_index": global_idx,
                        "position": snp_positions[global_idx],
                        "cn_change": f"{prev_cn}->{cn}",
                    }
            prev_cn = cn
    
    return None


# =============================================================================
# Utility Functions
# =============================================================================

def summarize_sv_result(sv_result: OtoaSvResult) -> Dict:
    """Convert SV result to dictionary for output."""
    return {
        "sv_type": sv_result.sv_type,
        "confidence": sv_result.confidence,
        "breakpoint": sv_result.breakpoint_region,
        "cn_pattern": sv_result.cn_pattern,
        "details": sv_result.details,
    }


def get_sv_description(sv_type: str) -> str:
    """Get human-readable description of SV type."""
    descriptions = {
        "normal": "Normal diploid (2 copies of OTOA, 2 copies of OTOAP1)",
        "cn2": "Normal diploid",
        "true_gene_del_het": "Heterozygous deletion of OTOA (1 copy remaining)",
        "true_gene_del_hom": "Homozygous deletion of OTOA (0 copies)",
        "true_gene_dup": "Duplication of OTOA (3 copies)",
        "true_gene_cn3": "OTOA copy number = 3",
        "true_gene_cn4": "OTOA copy number = 4",
        "pseudogene_del_het": "Heterozygous deletion of OTOAP1 pseudogene",
        "pseudogene_del_hom": "Homozygous deletion of OTOAP1 pseudogene",
        "partial_deletion": "Partial deletion affecting part of OTOA",
        "partial_duplication": "Partial duplication affecting part of OTOA",
        "gene_conversion": "Gene conversion between OTOA and OTOAP1",
        "combined_deletion": "Deletion affecting both OTOA and OTOAP1",
        "combined_duplication": "Duplication affecting both OTOA and OTOAP1",
        "hybrid": "Hybrid/fusion gene between OTOA and OTOAP1",
    }
    return descriptions.get(sv_type, f"Unknown SV type: {sv_type}")


def is_pathogenic_sv(sv_type: str) -> bool:
    """Determine if the SV is likely pathogenic for hearing loss."""
    pathogenic_types = [
        "true_gene_del_het",
        "true_gene_del_hom",
        "combined_deletion",
        "gene_conversion",
        "partial_deletion",
    ]
    return sv_type in pathogenic_types
