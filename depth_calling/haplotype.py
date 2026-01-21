#!/usr/bin/env python3
#
# Haplotype determination for OTOA variants
# Author: Yu-Hui Lin <yhlin.md05@nycu.edu.tw>
# OTOAlyzer - Haplotype phasing module
#
# Based on Cyrius (Illumina) and CyriPanel haplotype logic
# Original Cyrius Copyright (c) 2019-2020 Illumina, Inc.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""
Haplotype Determination Module for OTOAlyzer

This module implements read-based phasing to determine whether called variants
are on the same haplotype (cis) or different haplotypes (trans).

The approach follows Cyrius/CyriPanel methodology:
1. Extract bases from reads at target SNP positions
2. Translate bases to haplotype codes ("1" = true gene, "2" = pseudogene)
3. Group reads by haplotype pattern
4. Phase variants by checking which reads carry which variants

Key concepts:
- Haplotype string: e.g., "1121" means positions 0,1,3 have true gene allele,
  position 2 has pseudogene allele
- Cis variants: Both variants on same haplotype (same read/read-pair)
- Trans variants: Variants on different haplotypes (different reads)
"""

import logging
from collections import defaultdict
from typing import Dict, List, Optional, Set, Tuple, Any
import pysam

from .snp_count import passing_read


# =============================================================================
# Constants
# =============================================================================

MIN_MAPQ_HAPLOTYPE = 10  # Minimum mapping quality for haplotype analysis
MIN_READS_FOR_PHASING = 2  # Minimum reads required for confident phasing
MIN_DIFF_SITES_FOR_CLASSIFICATION = 2  # Minimum differentiating sites to classify a read
PSEUDOGENE_ALLELE_THRESHOLD = 0.7  # Fraction of pseudogene alleles to classify as pseudogene read


# =============================================================================
# Pseudogene Read Filtering
# =============================================================================

def classify_read_origin(
    read: pysam.AlignedSegment,
    snp_db: Any,
    min_sites: int = MIN_DIFF_SITES_FOR_CLASSIFICATION,
    pseudo_threshold: float = PSEUDOGENE_ALLELE_THRESHOLD,
) -> str:
    """
    Classify whether a read originates from true gene or pseudogene.

    Uses differentiating SNP sites to determine read origin. If a read
    carries mostly pseudogene-specific alleles, it's likely a mis-mapped
    pseudogene read.

    Args:
        read: pysam AlignedSegment object
        snp_db: SNP database from get_snp_position() containing:
            - dsnp1: True gene SNP dictionary {pos_idx: "true_base_pseudo_base"}
            - nchr: Chromosome name
        min_sites: Minimum differentiating sites required for classification
        pseudo_threshold: Fraction of pseudogene alleles to classify as pseudogene

    Returns:
        "true_gene", "pseudogene", or "unknown"
    """
    if snp_db is None:
        return "unknown"

    true_gene_count = 0
    pseudogene_count = 0

    # Get read sequence and reference positions
    read_seq = read.query_sequence
    if read_seq is None:
        return "unknown"

    # Get aligned pairs (query_pos, ref_pos)
    aligned_pairs = read.get_aligned_pairs(matches_only=True)
    if not aligned_pairs:
        return "unknown"

    # Build a lookup of reference position -> query position
    ref_to_query = {ref_pos: query_pos for query_pos, ref_pos in aligned_pairs}

    # Check each differentiating SNP site
    for snp_key, alleles in snp_db.dsnp1.items():
        # Extract position from key (format: "position_index")
        try:
            snp_pos = int(snp_key.split("_")[0])
        except (ValueError, IndexError):
            continue

        # Check if read covers this position (0-based in pysam)
        ref_pos_0based = snp_pos - 1
        if ref_pos_0based not in ref_to_query:
            continue

        query_pos = ref_to_query[ref_pos_0based]
        if query_pos >= len(read_seq):
            continue

        read_base = read_seq[query_pos].upper()

        # Parse alleles: "true_gene_base_pseudogene_base"
        try:
            true_allele, pseudo_allele = alleles.split("_")
        except ValueError:
            continue

        # Compare read base to expected alleles
        if read_base == true_allele.upper():
            true_gene_count += 1
        elif read_base == pseudo_allele.upper():
            pseudogene_count += 1

    total_sites = true_gene_count + pseudogene_count

    if total_sites < min_sites:
        return "unknown"

    pseudo_fraction = pseudogene_count / total_sites

    if pseudo_fraction >= pseudo_threshold:
        return "pseudogene"
    elif pseudo_fraction <= (1 - pseudo_threshold):
        return "true_gene"
    else:
        return "unknown"


def get_pseudogene_reads(
    bamfile_handle: pysam.AlignmentFile,
    snp_db: Any,
    region_chrom: str,
    region_start: int,
    region_end: int,
    min_mapq: int = MIN_MAPQ_HAPLOTYPE,
) -> Set[str]:
    """
    Identify reads that likely originate from the pseudogene.

    Scans all reads in the region and classifies them based on
    differentiating SNP sites.

    Args:
        bamfile_handle: Open pysam AlignmentFile
        snp_db: SNP database from get_snp_position()
        region_chrom: Chromosome name
        region_start: Region start position
        region_end: Region end position
        min_mapq: Minimum mapping quality

    Returns:
        Set of read names classified as pseudogene origin
    """
    pseudogene_read_names = set()

    for read in bamfile_handle.fetch(region_chrom, region_start, region_end):
        if read.mapping_quality < min_mapq:
            continue
        if read.is_secondary or read.is_supplementary or read.is_duplicate:
            continue

        origin = classify_read_origin(read, snp_db)
        if origin == "pseudogene":
            pseudogene_read_names.add(read.query_name)

    return pseudogene_read_names


# =============================================================================
# Core Haplotype Functions (adapted from Cyrius/CyriPanel)
# =============================================================================

def get_bases_per_read(
    bamfile_handle: pysam.AlignmentFile,
    base_db: Any,
    target_positions: List[int],
    region: Optional[int] = None,
    min_mapq: int = MIN_MAPQ_HAPLOTYPE,
) -> Dict[str, Dict[int, Optional[str]]]:
    """
    Extract bases from reads at target SNP positions.

    For each read that covers target positions, record the base observed
    at each position. This allows tracking which alleles co-occur on the
    same read (i.e., same haplotype).

    Args:
        bamfile_handle: Open pysam AlignmentFile
        base_db: SNP database from get_snp_position() containing:
            - dsnp1: True gene SNP dictionary
            - dsnp2: Pseudogene SNP dictionary
            - nchr: Chromosome name
            - dindex: Position index dictionary
        target_positions: List of SNP indices to analyze
        region: 0 for true gene only, 1 for pseudogene only, None for both
        min_mapq: Minimum mapping quality

    Returns:
        Dictionary: {read_name: {position_index: base_or_None}}
    """
    dread = {}
    nchr = base_db.nchr
    dindex = base_db.dindex

    # Select which SNP dictionaries to use based on region
    dsnps = [base_db.dsnp1, base_db.dsnp2]
    if region is not None:
        if region == 0:
            dsnps = [base_db.dsnp1]
        elif region == 1:
            dsnps = [base_db.dsnp2]

    for dsnp in dsnps:
        for snp_position_ori in dsnp:
            dsnp_index = dindex[snp_position_ori]
            snp_position = int(snp_position_ori.split("_")[0])

            # Only process target positions
            if dsnp_index not in target_positions:
                continue

            for pileupcolumn in bamfile_handle.pileup(
                nchr,
                snp_position - 1,
                snp_position + 1,
                truncate=True,
                stepper="nofilter",
                ignore_overlaps=False,
                ignore_orphan=False,
            ):
                site_position = pileupcolumn.pos + 1
                if site_position != snp_position:
                    continue

                reg1_allele, reg2_allele = dsnp[snp_position_ori].split("_")

                for read in pileupcolumn.pileups:
                    if not passing_read(read):
                        continue
                    if read.alignment.mapping_quality < min_mapq:
                        continue

                    read_name = read.alignment.query_name
                    read_seq = read.alignment.query_sequence
                    start_pos = read.query_position
                    end_pos = start_pos + min(len(reg1_allele), len(reg2_allele))

                    if end_pos >= len(read_seq):
                        continue

                    hap = read_seq[start_pos:end_pos]

                    # Initialize read entry if new
                    if read_name not in dread:
                        dread[read_name] = {}
                        for pos in target_positions:
                            dread[read_name][pos] = None

                    # Handle conflicting bases (mark as None)
                    if dread[read_name][dsnp_index] not in [None, hap]:
                        dread[read_name][dsnp_index] = None
                    else:
                        dread[read_name][dsnp_index] = hap

    return dread


def get_base1_base2(
    base_db: Any,
    target_positions: List[int],
) -> Tuple[List[str], List[str]]:
    """
    Get reference alleles for true gene (base1) and pseudogene (base2)
    at target positions.

    Args:
        base_db: SNP database from get_snp_position()
        target_positions: List of SNP indices

    Returns:
        Tuple of (base1_list, base2_list) for true gene and pseudogene alleles
    """
    base1 = []
    base2 = []

    # Sort positions to ensure consistent ordering
    sorted_positions = sorted(target_positions)

    for pos in sorted_positions:
        # Find the SNP entry for this position
        found = False
        for snp_key in base_db.dsnp1:
            if base_db.dindex.get(snp_key) == pos:
                alleles = base_db.dsnp1[snp_key].split("_")
                base1.append(alleles[0])
                base2.append(alleles[1])
                found = True
                break

        if not found:
            base1.append("N")
            base2.append("N")

    return base1, base2


def get_hap_counts(
    dread: Dict[str, Dict[int, Optional[str]]],
    base1: List[str],
    base2: List[str],
    target_positions: List[int],
) -> Dict[str, str]:
    """
    Translate observed bases into haplotype strings.

    For each read, convert the observed bases to a haplotype code:
    - "1" = matches true gene allele (base1)
    - "2" = matches pseudogene allele (base2)
    - "x" = missing or unknown

    Args:
        dread: Dictionary of read bases from get_bases_per_read()
        base1: True gene reference alleles
        base2: Pseudogene reference alleles
        target_positions: List of target position indices

    Returns:
        Dictionary: {read_name: haplotype_string}
    """
    dread_hap = {}
    sorted_positions = sorted(target_positions)

    for read_name, read_bases in dread.items():
        pos_list = ["x"] * len(sorted_positions)

        for i, pos in enumerate(sorted_positions):
            base = read_bases.get(pos)

            if base is None:
                continue

            # Check if base matches true gene allele
            for allele in base1[i].split(","):
                if base.upper() == allele.upper():
                    pos_list[i] = "1"
                    break

            # Check if base matches pseudogene allele
            if pos_list[i] == "x":
                for allele in base2[i].split(","):
                    if base.upper() == allele.upper():
                        pos_list[i] = "2"
                        break

        dread_hap[read_name] = "".join(pos_list)

    return dread_hap


def extract_hap(
    dhaplotype: Dict[str, str],
    positions_to_extract: List[int],
) -> Dict[str, List[int]]:
    """
    Extract haplotypes at specific positions from the full haplotype strings.

    Only includes reads with complete data (no 'x') at the requested positions.

    Args:
        dhaplotype: Dictionary of {read_name: haplotype_string}
        positions_to_extract: List of position indices to extract

    Returns:
        Dictionary: {haplotype_pattern: [1, 1, ...]} (list of read counts)
    """
    hap_count = {}

    for read_name, hap in dhaplotype.items():
        # Extract bases at requested positions
        try:
            hap_base = [hap[pos] for pos in positions_to_extract]
        except IndexError:
            continue

        # Only include complete haplotypes (no missing data)
        if "x" not in hap_base:
            hap_pattern = "".join(hap_base)
            if hap_pattern not in hap_count:
                hap_count[hap_pattern] = []
            hap_count[hap_pattern].append(1)

    return hap_count


def get_haplotypes_from_bam_single_region(
    bamfile_handle: pysam.AlignmentFile,
    base_db: Any,
    target_positions: List[int],
) -> Dict[str, str]:
    """
    Get haplotype assignments for reads in the true gene region.

    This is the main entry point for haplotype extraction.

    Args:
        bamfile_handle: Open pysam AlignmentFile
        base_db: SNP database from get_snp_position()
        target_positions: List of SNP indices to analyze

    Returns:
        Dictionary: {read_name: haplotype_string}
    """
    # Extract bases from reads at target positions (true gene region only)
    dread = get_bases_per_read(
        bamfile_handle, base_db, target_positions, region=0, min_mapq=MIN_MAPQ_HAPLOTYPE
    )

    # Get reference alleles
    base1, base2 = get_base1_base2(base_db, target_positions)

    # Convert to haplotype strings
    dhaplotype = get_hap_counts(dread, base1, base2, target_positions)

    return dhaplotype


# =============================================================================
# Variant Phasing Functions
# =============================================================================

def phase_variants(
    bamfile_handle: pysam.AlignmentFile,
    variants: List[Any],
    chromosome: str,
    snp_db: Optional[Any] = None,
    min_reads: int = MIN_READS_FOR_PHASING,
) -> Dict[str, Dict]:
    """
    Phase called variants to determine cis/trans relationships.

    Uses read-based phasing to determine if pairs of variants occur on
    the same haplotype (cis) or different haplotypes (trans).

    Args:
        bamfile_handle: Open pysam AlignmentFile
        variants: List of OtoaVariantCall objects
        chromosome: Chromosome name
        snp_db: SNP database for pseudogene read filtering (optional but recommended)
        min_reads: Minimum reads required for confident phasing

    Returns:
        Dictionary with phasing information:
        {
            "variant_haplotypes": {variant_id: haplotype_assignment},
            "phase_sets": [[var1, var2, ...], ...],  # Variants in same phase set
            "phase_relationships": {(var1, var2): "cis"|"trans"|"unknown"},
            "pseudogene_reads_filtered": int,  # Number of pseudogene reads excluded
        }
    """
    if len(variants) < 2:
        return {
            "variant_haplotypes": {},
            "phase_sets": [],
            "phase_relationships": {},
            "pseudogene_reads_filtered": 0,
        }

    # Get variant positions
    variant_positions = [(v.position, v) for v in variants]
    variant_positions.sort(key=lambda x: x[0])

    # Identify pseudogene reads in the variant region (if snp_db provided)
    pseudogene_read_names = set()
    if snp_db is not None:
        min_pos = min(v.position for v in variants)
        max_pos = max(v.position for v in variants)
        pseudogene_read_names = get_pseudogene_reads(
            bamfile_handle,
            snp_db,
            chromosome,
            min_pos - 1,
            max_pos + 1,
        )
        if pseudogene_read_names:
            logging.debug(f"Identified {len(pseudogene_read_names)} pseudogene reads to filter from phasing")

    # Extract reads covering variant positions (excluding pseudogene reads)
    reads_per_variant = {}
    variant_alleles_per_read = defaultdict(dict)
    filtered_count = 0

    for var in variants:
        reads_per_variant[var.variant_id] = {"ref": set(), "alt": set()}

        for pileupcolumn in bamfile_handle.pileup(
            chromosome,
            var.position - 1,
            var.position,
            truncate=True,
            min_mapping_quality=MIN_MAPQ_HAPLOTYPE,
        ):
            if pileupcolumn.pos + 1 != var.position:
                continue

            for pileup_read in pileupcolumn.pileups:
                if pileup_read.is_del or pileup_read.is_refskip:
                    continue

                try:
                    read_name = pileup_read.alignment.query_name

                    # Skip reads identified as pseudogene origin
                    if read_name in pseudogene_read_names:
                        filtered_count += 1
                        continue

                    base = pileup_read.alignment.query_sequence[pileup_read.query_position].upper()

                    if base == var.ref_base.upper():
                        reads_per_variant[var.variant_id]["ref"].add(read_name)
                        variant_alleles_per_read[read_name][var.variant_id] = "ref"
                    elif base == var.alt_base.upper():
                        reads_per_variant[var.variant_id]["alt"].add(read_name)
                        variant_alleles_per_read[read_name][var.variant_id] = "alt"
                except (TypeError, IndexError):
                    continue

    if filtered_count > 0:
        logging.info(f"Filtered {filtered_count} pseudogene read occurrences from phasing analysis")

    # Determine phase relationships between variant pairs
    phase_relationships = {}

    for i, var1 in enumerate(variants):
        for var2 in variants[i+1:]:
            relationship = determine_phase_relationship(
                var1, var2, variant_alleles_per_read, min_reads
            )
            phase_relationships[(var1.variant_id, var2.variant_id)] = relationship

    # Build phase sets (groups of variants on same haplotype)
    phase_sets = build_phase_sets(variants, phase_relationships)

    # Assign haplotype labels
    variant_haplotypes = assign_haplotype_labels(variants, phase_sets)

    return {
        "variant_haplotypes": variant_haplotypes,
        "phase_sets": phase_sets,
        "phase_relationships": phase_relationships,
        "pseudogene_reads_filtered": len(pseudogene_read_names),
    }


def determine_phase_relationship(
    var1: Any,
    var2: Any,
    variant_alleles_per_read: Dict[str, Dict[str, str]],
    min_reads: int,
) -> str:
    """
    Determine if two variants are in cis (same haplotype) or trans (different haplotypes).

    Logic:
    - If reads carry ALT for both variants → cis
    - If reads carry ALT for one and REF for other → trans
    - If mixed or insufficient reads → unknown

    Args:
        var1, var2: Variant objects
        variant_alleles_per_read: {read_name: {variant_id: "ref"|"alt"}}
        min_reads: Minimum reads for confident call

    Returns:
        "cis", "trans", or "unknown"
    """
    cis_count = 0  # Both ALT or both REF on same read
    trans_count = 0  # ALT/REF mismatch on same read

    for read_name, alleles in variant_alleles_per_read.items():
        if var1.variant_id not in alleles or var2.variant_id not in alleles:
            continue

        allele1 = alleles[var1.variant_id]
        allele2 = alleles[var2.variant_id]

        if allele1 == allele2:
            # Both ref or both alt = cis configuration
            cis_count += 1
        else:
            # One ref, one alt = trans configuration
            trans_count += 1

    total = cis_count + trans_count

    if total < min_reads:
        return "unknown"

    # Require clear majority for confident call
    if cis_count > trans_count * 2:
        return "cis"
    elif trans_count > cis_count * 2:
        return "trans"
    else:
        return "unknown"


def build_phase_sets(
    variants: List[Any],
    phase_relationships: Dict[Tuple[str, str], str],
) -> List[List[str]]:
    """
    Build phase sets (groups of variants on the same haplotype).

    Uses union-find to group variants that are in cis.

    Args:
        variants: List of variant objects
        phase_relationships: {(var1_id, var2_id): "cis"|"trans"|"unknown"}

    Returns:
        List of phase sets, each containing variant IDs
    """
    # Initialize each variant in its own set
    parent = {v.variant_id: v.variant_id for v in variants}

    def find(x):
        if parent[x] != x:
            parent[x] = find(parent[x])
        return parent[x]

    def union(x, y):
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py

    # Union variants that are in cis
    for (var1_id, var2_id), relationship in phase_relationships.items():
        if relationship == "cis":
            union(var1_id, var2_id)

    # Build phase sets
    sets = defaultdict(list)
    for v in variants:
        sets[find(v.variant_id)].append(v.variant_id)

    return list(sets.values())


def assign_haplotype_labels(
    variants: List[Any],
    phase_sets: List[List[str]],
) -> Dict[str, str]:
    """
    Assign haplotype labels (H1, H2, etc.) to variants based on phase sets.

    Args:
        variants: List of variant objects
        phase_sets: List of phase sets from build_phase_sets()

    Returns:
        {variant_id: "H1"|"H2"|...}
    """
    haplotypes = {}

    for i, phase_set in enumerate(phase_sets):
        hap_label = f"H{i + 1}"
        for var_id in phase_set:
            haplotypes[var_id] = hap_label

    return haplotypes


# =============================================================================
# Compound Heterozygosity Detection
# =============================================================================

def detect_compound_het(
    phase_relationships: Dict[Tuple[str, str], str],
    variants: List[Any],
) -> List[Tuple[str, str]]:
    """
    Detect compound heterozygous variant pairs (trans configuration).

    Compound heterozygotes have one variant on each haplotype, which can be
    clinically significant when both affect the same gene.

    Args:
        phase_relationships: {(var1_id, var2_id): "cis"|"trans"|"unknown"}
        variants: List of variant objects

    Returns:
        List of (var1_id, var2_id) tuples that are compound het
    """
    compound_het_pairs = []

    for (var1_id, var2_id), relationship in phase_relationships.items():
        if relationship == "trans":
            compound_het_pairs.append((var1_id, var2_id))

    return compound_het_pairs


def summarize_phasing(
    phasing_result: Dict,
    variants: List[Any],
) -> Dict:
    """
    Generate summary statistics for phasing results.

    Args:
        phasing_result: Output from phase_variants()
        variants: List of variant objects

    Returns:
        Summary dictionary
    """
    relationships = phasing_result.get("phase_relationships", {})

    cis_count = sum(1 for r in relationships.values() if r == "cis")
    trans_count = sum(1 for r in relationships.values() if r == "trans")
    unknown_count = sum(1 for r in relationships.values() if r == "unknown")

    compound_het = detect_compound_het(relationships, variants)

    return {
        "total_variants": len(variants),
        "total_pairs": len(relationships),
        "cis_pairs": cis_count,
        "trans_pairs": trans_count,
        "unknown_pairs": unknown_count,
        "phase_sets": len(phasing_result.get("phase_sets", [])),
        "compound_het_pairs": len(compound_het),
        "compound_het_variants": compound_het,
        "pseudogene_reads_filtered": phasing_result.get("pseudogene_reads_filtered", 0),
    }
