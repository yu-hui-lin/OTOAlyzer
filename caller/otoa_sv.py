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

SNP Distribution in OTOA_SNP_38.txt (85 SNPs total):
- upstreamparalog: 8 SNPs
- upstreamparalog_150bp: 2 SNPs
- intron21: 46 SNPs (spans 5prime and middle regions)
- exon22: 2 SNPs
- intron22: 5 SNPs
- intron23: 1 SNP
- intron24: 1 SNP
- intron27: 5 SNPs
- intron28: 13 SNPs
- exon29: 1 SNP
- down_exon29: 1 SNP
Total: 85 SNPs

Regional Division (4 regions, defined by 'region' column in SNP file):
1. upstream:  upstreamparalog + upstreamparalog_150bp (10 SNPs, indices 0-9)
2. 5prime:    intron21 (5' portion) (33 SNPs, indices 10-42)
3. middle:    intron21 (3' portion) + exon22 + intron22 + intron23 (21 SNPs, indices 43-63)
4. 3prime:    intron24 + intron27 + intron28 + exon29 + down_exon29 (21 SNPs, indices 64-84)

SV Types Detected:
- Deletions (heterozygous/homozygous) of true gene
- Deletions of pseudogene
- Duplications of true gene
- Duplications of pseudogene
- Gene conversions (partial CN changes)
- Hybrid genes
"""

import logging
from collections import namedtuple, Counter, OrderedDict
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
        "upstream",      # upstream paralog region
        "region_5prime", # 5' portion of homologous region
        "region_middle", # middle portion of homologous region
        "region_3prime", # 3' portion of homologous region
    ]
)

# =============================================================================
# Region Boundary Computation
# =============================================================================

# Expected region labels in order (must match 'region' column in SNP file)
EXPECTED_REGION_ORDER = ["upstream", "5prime", "middle", "3prime"]

# Looser threshold when primary consensus fails
CONSENSUS_MIN_LOOSE = 5

# Minimum fraction of region sites for reliable consensus
CONSENSUS_MIN_FRACTION = 0.5


def compute_region_boundaries(dregion: Dict[int, str], num_sites: int) -> Dict[str, Tuple[int, int]]:
    """
    Compute region boundaries (start, end indices) from the region assignments
    parsed from the SNP file's 'region' column.

    This replaces hardcoded constants, making the code robust to changes in
    the number of differentiating sites or region borders.

    Args:
        dregion: Dict mapping SNP index to region label
                 (e.g., {0: "upstream", 1: "upstream", ..., 10: "5prime", ...})
        num_sites: Total number of SNP sites

    Returns:
        OrderedDict mapping region name to (start_index, end_index) tuple
        where end_index is exclusive (Python slice convention).
        e.g., {"upstream": (0, 10), "5prime": (10, 43), "middle": (43, 64), "3prime": (64, 85)}
    """
    if not dregion:
        logging.warning("No region information available; using equal-split fallback")
        quarter = num_sites // 4
        return OrderedDict([
            ("upstream", (0, quarter)),
            ("5prime", (quarter, 2 * quarter)),
            ("middle", (2 * quarter, 3 * quarter)),
            ("3prime", (3 * quarter, num_sites)),
        ])

    # Collect all indices per region label
    region_indices = {}
    for idx in sorted(dregion.keys()):
        label = dregion[idx]
        region_indices.setdefault(label, [])
        region_indices[label].append(idx)

    # Build boundaries: (min_index, max_index + 1) for each region in order
    boundaries = OrderedDict()
    for region_name in EXPECTED_REGION_ORDER:
        if region_name in region_indices:
            indices = region_indices[region_name]
            boundaries[region_name] = (min(indices), max(indices) + 1)
        else:
            logging.warning(f"Region '{region_name}' not found in SNP file")

    return boundaries


def get_min_sites_for_region(region_size: int, fraction: float = CONSENSUS_MIN_FRACTION) -> int:
    """
    Calculate minimum number of sites required for consensus in a region.

    Args:
        region_size: Number of SNP sites in the region
        fraction: Minimum fraction of sites required (default 0.5)

    Returns:
        Minimum site count (at least 1)
    """
    return max(1, int(region_size * fraction))


# =============================================================================
# Main SV Detection
# =============================================================================

def get_otoa_sv_call(
    total_cn: int,
    cn_call_per_site: List[Optional[int]],
    true_gene_cn: int,
    pseudogene_cn: int,
    region_boundaries: Optional[Dict[str, Tuple[int, int]]] = None,
) -> OtoaSvResult:
    """
    Determine SV type based on regional CN patterns.
    
    This is the main entry point for SV detection. It analyzes
    CN patterns across the OTOA gene regions to detect SVs.
    
    Args:
        total_cn: Total copy number (true gene + pseudogene)
        cn_call_per_site: Per-site CN calls from Poisson model
        true_gene_cn: CN of true OTOA gene (from OTOA_unique)
        pseudogene_cn: CN of pseudogene (total_cn - true_gene_cn)
        region_boundaries: Dict mapping region names to (start, end) index tuples.
                           Computed dynamically from SNP file if provided.
        
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
    consensus = get_regional_consensus(cn_call_per_site, total_cn, region_boundaries)
    
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
    region_boundaries: Optional[Dict[str, Tuple[int, int]]] = None,
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
        region_boundaries: Dict mapping region names to (start, end) index tuples.
                           If None, falls back to equal-split partitioning.
        
    Returns:
        OtoaRegionalCN named tuple with consensus CN for each region
    """
    # Handle empty or short input
    if not cn_call_per_site:
        return OtoaRegionalCN(None, None, None, None)
    
    num_sites = len(cn_call_per_site)

    # Use provided boundaries or compute fallback
    if region_boundaries is None:
        logging.warning("No region boundaries provided; using equal-split fallback")
        quarter = num_sites // 4
        region_boundaries = OrderedDict([
            ("upstream", (0, quarter)),
            ("5prime", (quarter, 2 * quarter)),
            ("middle", (2 * quarter, 3 * quarter)),
            ("3prime", (3 * quarter, num_sites)),
        ])

    def _extract_region_sites(start, end):
        """Extract non-None CN values for a region."""
        return [a for a in cn_call_per_site[start:end] if a is not None]

    # Extract sites per region
    region_sites = {}
    for region_name, (start, end) in region_boundaries.items():
        region_sites[region_name] = _extract_region_sites(start, end)
        region_size = end - start
        logging.debug(f"Region {region_name}: {len(region_sites[region_name])}/{region_size} valid sites "
                      f"(indices {start}-{end-1})")

    # Calculate consensus for each region
    upstream_consensus = _get_region_consensus(
        region_sites.get("upstream", []),
        get_min_sites_for_region(
            region_boundaries["upstream"][1] - region_boundaries["upstream"][0]
        ) if "upstream" in region_boundaries else 1,
        total_cn,
    )

    consensus_5prime = _get_region_consensus(
        region_sites.get("5prime", []),
        get_min_sites_for_region(
            region_boundaries["5prime"][1] - region_boundaries["5prime"][0]
        ) if "5prime" in region_boundaries else 1,
        total_cn,
    )

    consensus_middle = _get_region_consensus(
        region_sites.get("middle", []),
        get_min_sites_for_region(
            region_boundaries["middle"][1] - region_boundaries["middle"][0]
        ) if "middle" in region_boundaries else 1,
        total_cn,
    )

    consensus_3prime = _get_region_consensus(
        region_sites.get("3prime", []),
        get_min_sites_for_region(
            region_boundaries["3prime"][1] - region_boundaries["3prime"][0]
        ) if "3prime" in region_boundaries else 1,
        total_cn,
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
        true_gene_cn: True gene CN (from OTOA_unique)
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
    
    # --- Simple whole-gene events for TRUE GENE (OTOA) ---
    
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
    
    # --- Simple whole-gene events for PSEUDOGENE (OTOAP1) ---
    
    # Heterozygous deletion of pseudogene (CN=1)
    if true_gene_cn == 2 and pseudogene_cn == 1:
        return ("pseudogene_del_het", 0.85, ["pseudogene_deletion"])
    
    # Homozygous deletion of pseudogene (CN=0)
    if true_gene_cn == 2 and pseudogene_cn == 0:
        return ("pseudogene_del_hom", 0.90, ["pseudogene_deletion_hom"])
    
    # Duplication of pseudogene (CN=3)
    if true_gene_cn == 2 and pseudogene_cn == 3:
        return ("pseudogene_dup", 0.85, ["pseudogene_duplication"])
    
    # Triplication or higher of pseudogene (CN>=4)
    if true_gene_cn == 2 and pseudogene_cn >= 4:
        return (f"pseudogene_cn{pseudogene_cn}", 0.80, ["pseudogene_multiplication"])
    
    # --- Combined events: both true gene and pseudogene affected ---
    
    # True gene deletion + pseudogene duplication
    if true_gene_cn == 1 and pseudogene_cn == 3:
        return ("true_del_pseudo_dup", 0.75, ["true_gene_deletion", "pseudogene_duplication"])
    
    if true_gene_cn == 1 and pseudogene_cn >= 4:
        return (f"true_del_pseudo_cn{pseudogene_cn}", 0.70, ["true_gene_deletion", "pseudogene_multiplication"])
    
    # True gene duplication + pseudogene deletion
    if true_gene_cn == 3 and pseudogene_cn == 1:
        return ("true_dup_pseudo_del", 0.75, ["true_gene_duplication", "pseudogene_deletion"])
    
    if true_gene_cn >= 4 and pseudogene_cn == 1:
        return (f"true_cn{true_gene_cn}_pseudo_del", 0.70, ["true_gene_multiplication", "pseudogene_deletion"])
    
    # True gene duplication + pseudogene duplication (both increased)
    if true_gene_cn == 3 and pseudogene_cn == 3:
        return ("both_dup", 0.75, ["true_gene_duplication", "pseudogene_duplication"])
    
    if true_gene_cn >= 3 and pseudogene_cn >= 3:
        return (f"true_cn{true_gene_cn}_pseudo_cn{pseudogene_cn}", 0.70, 
                ["true_gene_multiplication", "pseudogene_multiplication"])
    
    # True gene deletion + pseudogene deletion (both decreased)
    if true_gene_cn == 1 and pseudogene_cn == 1:
        return ("both_del_het", 0.75, ["true_gene_deletion", "pseudogene_deletion"])
    
    if true_gene_cn == 0 and pseudogene_cn == 1:
        return ("true_del_hom_pseudo_del_het", 0.80, ["true_gene_deletion_hom", "pseudogene_deletion"])
    
    if true_gene_cn == 1 and pseudogene_cn == 0:
        return ("true_del_het_pseudo_del_hom", 0.80, ["true_gene_deletion", "pseudogene_deletion_hom"])
    
    if true_gene_cn == 0 and pseudogene_cn == 0:
        return ("both_del_hom", 0.85, ["true_gene_deletion_hom", "pseudogene_deletion_hom"])
    
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
    
    # Both true gene and pseudogene affected (not caught above)
    if true_gene_cn != 2 and pseudogene_cn != 2:
        if true_gene_cn < 2 and pseudogene_cn < 2:
            return ("combined_deletion", 0.70, ["both_deleted"])
        elif true_gene_cn > 2 and pseudogene_cn > 2:
            return ("combined_duplication", 0.65, ["both_duplicated"])
        else:
            # One increased, one decreased - complex rearrangement
            return ("hybrid", 0.60, ["complex_rearrangement"])
    
    # Check for gene conversion when CNV caller shows normal but regions vary
    if true_gene_cn == 2 and pseudogene_cn == 2:
        if has_regional_variation(consensus):
            return ("gene_conversion", 0.70, ["regional_variation"])
    
    # Cannot determine - provide informative SV type based on available data
    if true_gene_cn != 2 or pseudogene_cn != 2:
        return (f"unclassified_true{true_gene_cn}_pseudo{pseudogene_cn}", 0.50, 
                ["unclassified_cn_change"])
    
    logging.warning(f"Could not classify SV: true_cn={true_gene_cn}, pseudo_cn={pseudogene_cn}, "
                    f"consensus={consensus}")
    return ("normal", 0.90, [])


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
        "del_upstream_to_5prime": "upstream_intron21_5prime_boundary",
        "dup_upstream_to_5prime": "upstream_intron21_5prime_boundary",
        "del_5prime_to_middle": "intron21_middle_boundary",
        "dup_5prime_to_middle": "intron21_middle_boundary",
        "del_middle_to_3prime": "intron23_intron24_boundary",
        "dup_middle_to_3prime": "intron23_intron24_boundary",
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
        # Normal
        "normal": "Normal diploid (2 copies of OTOA, 2 copies of OTOAP1)",
        "cn2": "Normal diploid",
        
        # True gene (OTOA) events
        "true_gene_del_het": "Heterozygous deletion of OTOA (1 copy remaining)",
        "true_gene_del_hom": "Homozygous deletion of OTOA (0 copies)",
        "true_gene_dup": "Duplication of OTOA (3 copies)",
        "true_gene_cn3": "OTOA copy number = 3",
        "true_gene_cn4": "OTOA copy number = 4",
        "true_gene_cn5": "OTOA copy number = 5",
        "true_gene_cn6": "OTOA copy number = 6",
        
        # Pseudogene (OTOAP1) events
        "pseudogene_del_het": "Heterozygous deletion of OTOAP1 pseudogene (1 copy remaining)",
        "pseudogene_del_hom": "Homozygous deletion of OTOAP1 pseudogene (0 copies)",
        "pseudogene_dup": "Duplication of OTOAP1 pseudogene (3 copies)",
        "pseudogene_cn3": "OTOAP1 pseudogene copy number = 3",
        "pseudogene_cn4": "OTOAP1 pseudogene copy number = 4",
        "pseudogene_cn5": "OTOAP1 pseudogene copy number = 5",
        "pseudogene_cn6": "OTOAP1 pseudogene copy number = 6",
        
        # Combined events
        "true_del_pseudo_dup": "OTOA deletion + OTOAP1 duplication",
        "true_dup_pseudo_del": "OTOA duplication + OTOAP1 deletion",
        "both_dup": "Duplication of both OTOA and OTOAP1",
        "both_del_het": "Heterozygous deletion of both OTOA and OTOAP1",
        "both_del_hom": "Homozygous deletion of both OTOA and OTOAP1",
        "true_del_hom_pseudo_del_het": "Homozygous OTOA deletion + heterozygous OTOAP1 deletion",
        "true_del_het_pseudo_del_hom": "Heterozygous OTOA deletion + homozygous OTOAP1 deletion",
        
        # Partial/complex events
        "partial_deletion": "Partial deletion affecting part of OTOA",
        "partial_duplication": "Partial duplication affecting part of OTOA",
        "gene_conversion": "Gene conversion between OTOA and OTOAP1",
        "combined_deletion": "Deletion affecting both OTOA and OTOAP1",
        "combined_duplication": "Duplication affecting both OTOA and OTOAP1",
        "hybrid": "Hybrid/fusion gene between OTOA and OTOAP1",
    }
    
    # Check for exact match first
    if sv_type in descriptions:
        return descriptions[sv_type]
    
    # Handle dynamic CN types (e.g., "true_gene_cn7", "pseudogene_cn8")
    if sv_type and sv_type.startswith("true_gene_cn"):
        try:
            cn = sv_type.replace("true_gene_cn", "")
            return f"OTOA copy number = {cn}"
        except:
            pass
    
    if sv_type and sv_type.startswith("pseudogene_cn"):
        try:
            cn = sv_type.replace("pseudogene_cn", "")
            return f"OTOAP1 pseudogene copy number = {cn}"
        except:
            pass
    
    # Handle combined dynamic types
    if sv_type and "true" in sv_type and "pseudo" in sv_type:
        return f"Complex CN change: {sv_type}"
    
    return f"Unknown SV type: {sv_type}"


def is_pathogenic_sv(sv_type: str) -> bool:
    """
    Determine if the SV is likely pathogenic for hearing loss.
    
    OTOA-related hearing loss (DFNB22) is autosomal recessive,
    so pathogenic SVs are those that reduce functional OTOA copies.
    
    Note: Pseudogene duplications are generally NOT pathogenic
    as they don't affect the functional gene.
    """
    pathogenic_types = [
        # True gene deletions (reduce functional copies)
        "true_gene_del_het",
        "true_gene_del_hom",
        
        # Combined deletions
        "combined_deletion",
        "both_del_het",
        "both_del_hom",
        "true_del_hom_pseudo_del_het",
        "true_del_het_pseudo_del_hom",
        
        # Gene conversion (may disrupt function)
        "gene_conversion",
        
        # Partial events affecting true gene
        "partial_deletion",
        
        # Hybrid genes (may be non-functional)
        "hybrid",
    ]
    
    if sv_type in pathogenic_types:
        return True
    
    # Also check for dynamic types with true gene deletion
    if sv_type and "true_del" in sv_type:
        return True
    
    return False


def is_benign_sv(sv_type: str) -> bool:
    """
    Determine if the SV is likely benign.
    
    Pseudogene-only events are generally benign as they don't
    affect the functional OTOA gene.
    """
    benign_types = [
        "normal",
        "cn2",
        "pseudogene_del_het",
        "pseudogene_del_hom",
        "pseudogene_dup",
    ]
    
    if sv_type in benign_types:
        return True
    
    # Pseudogene CN changes are benign
    if sv_type and sv_type.startswith("pseudogene_cn"):
        return True
    
    return False
