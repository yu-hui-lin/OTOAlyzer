#!/usr/bin/env python3
#
# OTOAlyzer: CNV/SV and variant caller for OTOA gene and pseudogene
# Module: otoa_variants.py - De novo SNP and indel variant calling
#
# This module performs de novo variant discovery across the OTOA true gene region,
# NOT relying on any predefined variant database.

"""
De novo variant calling in OTOA true gene region.

Key features:
1. Scans entire OTOA true gene region for variants
2. Does NOT require predefined variant positions
3. Filters out pseudogene contamination using differentiating SNP sites
4. Reports all detected variants with quality metrics

The challenge: OTOA has a highly homologous pseudogene, so we must carefully
distinguish true variants from:
- Pseudogene reads mismapping to true gene
- Sequencing/alignment artifacts
- True somatic/germline variants
"""

import logging
from collections import namedtuple, defaultdict
from typing import List, Optional, Tuple, Dict, Set
from scipy.stats import fisher_exact
import math
import pysam

# =============================================================================
# Data Structures
# =============================================================================

OtoaVariantCall = namedtuple(
    "OtoaVariantCall",
    [
        "variant_id",       # Unique identifier (chr:pos:ref>alt)
        "chromosome",       # Chromosome
        "position",         # Genomic position (1-based)
        "ref_base",         # Reference base(s)
        "alt_base",         # Alternate base(s)
        "variant_type",     # SNV, insertion, deletion
        "allele_fraction",  # Variant allele fraction
        "total_depth",      # Total read depth at position
        "ref_reads",        # Reference allele reads
        "alt_reads",        # Alternate allele reads
        "zygosity",         # HET, HOM, or HEMI
        "filter_status",    # PASS or filter reason(s)
        "quality_score",    # Phred-scaled quality score
        "strand_bias_pval", # Strand bias p-value
        "in_homology_region", # Whether position is in high-homology region
        "pseudogene_evidence", # Evidence of pseudogene contamination
        "haplotype",        # Haplotype assignment (H1, H2, or None)
        "phase_set",        # Phase set ID for linked variants
    ]
)

# =============================================================================
# Constants and Thresholds
# =============================================================================

# Minimum thresholds for variant calling
MIN_TOTAL_DEPTH = 10        # Minimum total reads at position
MIN_ALT_READS = 3           # Minimum alternate allele reads
MIN_ALT_FRACTION = 0.10     # Minimum allele fraction to consider
MIN_BASE_QUALITY = 20       # Minimum base quality score
MIN_MAPPING_QUALITY = 20    # Minimum mapping quality

# Zygosity thresholds (for diploid regions)
HET_AF_MIN = 0.20           # Minimum AF for heterozygous
HET_AF_MAX = 0.80           # Maximum AF for heterozygous  
HOM_AF_MIN = 0.80           # Minimum AF for homozygous

# Quality filters
STRAND_BIAS_PVALUE = 0.001  # Strand bias significance threshold
MIN_STRAND_RATIO = 0.10     # Minimum reads on minor strand

# Pseudogene contamination detection
PSEUDOGENE_CONTAM_THRESHOLD = 0.15  # Max fraction suggesting contamination


# =============================================================================
# Main De Novo Variant Calling
# =============================================================================

def call_variants_denovo(
    bamfile: pysam.AlignmentFile,
    true_gene_region: Dict[str, any],
    pseudogene_region: Optional[Dict[str, any]],
    homology_sites: Optional[Set[int]],
    true_gene_cn: int,
    reference_fasta: Optional[str] = None,
) -> List[OtoaVariantCall]:
    """
    Perform de novo variant calling across the entire OTOA true gene region.
    
    This function scans every position in the true gene region and identifies
    variants without relying on a predefined variant database.
    
    Args:
        bamfile: Open pysam AlignmentFile object
        true_gene_region: Dict with 'chrom', 'start', 'end' for true gene
        pseudogene_region: Dict with 'chrom', 'start', 'end' for pseudogene (optional)
        homology_sites: Set of positions known to be in high-homology regions
        true_gene_cn: Copy number of true gene (affects zygosity calling)
        reference_fasta: Path to reference FASTA (optional, for indel calling)
        
    Returns:
        List of OtoaVariantCall objects for all detected variants
    """
    variants = []
    
    chrom = true_gene_region['chrom']
    start = true_gene_region['start']
    end = true_gene_region['end']
    
    logging.info(f"Scanning {chrom}:{start}-{end} for de novo variants...")
    
    # Get reference sequence if available
    ref_sequence = None
    if reference_fasta:
        try:
            ref_fasta = pysam.FastaFile(reference_fasta)
            ref_sequence = ref_fasta.fetch(chrom, start, end).upper()
            ref_fasta.close()
        except Exception as e:
            logging.warning(f"Could not load reference: {e}")
    
    # Track statistics
    positions_scanned = 0
    variants_found = 0
    
    # Iterate through pileup
    for pileup_column in bamfile.pileup(
        chrom,
        start,
        end,
        truncate=True,
        min_mapping_quality=MIN_MAPPING_QUALITY,
        min_base_quality=MIN_BASE_QUALITY,
        stepper='samtools',
    ):
        positions_scanned += 1
        position = pileup_column.pos + 1  # Convert to 1-based
        
        # Get base counts at this position
        base_counts, strand_counts, indel_counts = count_bases_at_position(pileup_column)
        
        total_depth = sum(base_counts.values())
        
        # Skip low-depth positions
        if total_depth < MIN_TOTAL_DEPTH:
            continue
        
        # Determine reference base
        ref_base = None
        if ref_sequence:
            ref_idx = pileup_column.pos - start
            if 0 <= ref_idx < len(ref_sequence):
                ref_base = ref_sequence[ref_idx]
        
        # If no reference, infer from most common base (with caution)
        if ref_base is None:
            ref_base = max(base_counts.keys(), key=lambda x: base_counts[x]) if base_counts else 'N'
        
        # Check for SNVs
        for alt_base, alt_count in base_counts.items():
            if alt_base == ref_base:
                continue
            
            # Calculate allele fraction
            af = alt_count / total_depth
            
            # Apply minimum thresholds
            if alt_count < MIN_ALT_READS or af < MIN_ALT_FRACTION:
                continue
            
            # This is a potential variant - perform detailed analysis
            var_call = analyze_variant(
                chrom=chrom,
                position=position,
                ref_base=ref_base,
                alt_base=alt_base,
                alt_count=alt_count,
                total_depth=total_depth,
                ref_count=base_counts.get(ref_base, 0),
                strand_counts=strand_counts,
                true_gene_cn=true_gene_cn,
                homology_sites=homology_sites,
                pseudogene_region=pseudogene_region,
                bamfile=bamfile,
            )
            
            if var_call is not None:
                variants.append(var_call)
                variants_found += 1
        
        # Check for indels
        for indel_key, indel_count in indel_counts.items():
            indel_type, indel_seq = indel_key
            
            af = indel_count / total_depth
            
            if indel_count < MIN_ALT_READS or af < MIN_ALT_FRACTION:
                continue
            
            # Analyze indel
            var_call = analyze_indel(
                chrom=chrom,
                position=position,
                ref_base=ref_base,
                indel_type=indel_type,
                indel_seq=indel_seq,
                indel_count=indel_count,
                total_depth=total_depth,
                true_gene_cn=true_gene_cn,
                homology_sites=homology_sites,
            )
            
            if var_call is not None:
                variants.append(var_call)
                variants_found += 1
    
    logging.info(f"Scanned {positions_scanned} positions, found {variants_found} variants")
    
    return variants


def count_bases_at_position(pileup_column) -> Tuple[Dict[str, int], Dict[str, Dict[str, int]], Dict[Tuple, int]]:
    """
    Count bases, strand information, and indels at a pileup position.
    
    Args:
        pileup_column: pysam PileupColumn object
        
    Returns:
        Tuple of (base_counts, strand_counts, indel_counts)
        - base_counts: {base: count}
        - strand_counts: {base: {'forward': n, 'reverse': n}}
        - indel_counts: {(type, seq): count}
    """
    base_counts = defaultdict(int)
    strand_counts = defaultdict(lambda: {'forward': 0, 'reverse': 0})
    indel_counts = defaultdict(int)
    
    for pileup_read in pileup_column.pileups:
        # Handle deletions
        if pileup_read.is_del:
            indel_counts[('del', str(pileup_read.indel))] += 1
            continue
        
        # Handle reference skips (splicing)
        if pileup_read.is_refskip:
            continue
        
        # Get the base
        try:
            base = pileup_read.alignment.query_sequence[pileup_read.query_position].upper()
        except (TypeError, IndexError):
            continue
        
        base_counts[base] += 1
        
        # Track strand
        if pileup_read.alignment.is_reverse:
            strand_counts[base]['reverse'] += 1
        else:
            strand_counts[base]['forward'] += 1
        
        # Check for insertions
        if pileup_read.indel > 0:
            # Get inserted sequence
            try:
                insert_start = pileup_read.query_position + 1
                insert_end = insert_start + pileup_read.indel
                insert_seq = pileup_read.alignment.query_sequence[insert_start:insert_end]
                indel_counts[('ins', insert_seq.upper())] += 1
            except (TypeError, IndexError):
                pass
    
    return dict(base_counts), dict(strand_counts), dict(indel_counts)


def analyze_variant(
    chrom: str,
    position: int,
    ref_base: str,
    alt_base: str,
    alt_count: int,
    total_depth: int,
    ref_count: int,
    strand_counts: Dict,
    true_gene_cn: int,
    homology_sites: Optional[Set[int]],
    pseudogene_region: Optional[Dict],
    bamfile: pysam.AlignmentFile,
) -> Optional[OtoaVariantCall]:
    """
    Perform detailed analysis of a potential SNV.
    
    Args:
        Various parameters describing the variant and context
        
    Returns:
        OtoaVariantCall if variant passes analysis, None otherwise
    """
    # Calculate allele fraction
    af = alt_count / total_depth
    
    # Calculate strand bias
    alt_strand = strand_counts.get(alt_base, {'forward': 0, 'reverse': 0})
    forward_alt = alt_strand['forward']
    reverse_alt = alt_strand['reverse']
    
    strand_bias_pval = calculate_strand_bias(forward_alt, reverse_alt, ref_count)
    
    # Check for pseudogene contamination
    pseudogene_evidence = None
    if pseudogene_region:
        pseudogene_evidence = check_pseudogene_signal(
            bamfile, chrom, position, alt_base, pseudogene_region
        )
    
    # Determine if in homology region
    in_homology = position in homology_sites if homology_sites else False
    
    # Determine zygosity
    zygosity = determine_zygosity(af, true_gene_cn)
    
    # Calculate quality score (Phred-scaled)
    quality_score = calculate_variant_quality(alt_count, total_depth, af)
    
    # Apply filters
    filter_status = apply_variant_filters(
        af=af,
        alt_count=alt_count,
        total_depth=total_depth,
        strand_bias_pval=strand_bias_pval,
        forward_alt=forward_alt,
        reverse_alt=reverse_alt,
        in_homology=in_homology,
        pseudogene_evidence=pseudogene_evidence,
    )
    
    # Create variant ID
    variant_id = f"{chrom}:{position}:{ref_base}>{alt_base}"
    
    return OtoaVariantCall(
        variant_id=variant_id,
        chromosome=chrom,
        position=position,
        ref_base=ref_base,
        alt_base=alt_base,
        variant_type="SNV",
        allele_fraction=round(af, 4),
        total_depth=total_depth,
        ref_reads=ref_count,
        alt_reads=alt_count,
        zygosity=zygosity,
        filter_status=filter_status,
        quality_score=round(quality_score, 1),
        strand_bias_pval=round(strand_bias_pval, 6) if strand_bias_pval else None,
        in_homology_region=in_homology,
        pseudogene_evidence=pseudogene_evidence,
        haplotype=None,  # Set later by phase_variants()
        phase_set=None,  # Set later by phase_variants()
    )


def analyze_indel(
    chrom: str,
    position: int,
    ref_base: str,
    indel_type: str,
    indel_seq: str,
    indel_count: int,
    total_depth: int,
    true_gene_cn: int,
    homology_sites: Optional[Set[int]],
) -> Optional[OtoaVariantCall]:
    """
    Analyze a potential insertion or deletion.
    
    Args:
        Various parameters describing the indel
        
    Returns:
        OtoaVariantCall if indel passes analysis, None otherwise
    """
    af = indel_count / total_depth
    
    # Determine ref/alt based on indel type
    if indel_type == 'ins':
        ref_display = ref_base
        alt_display = ref_base + indel_seq
        var_type = "INS"
    else:  # deletion
        ref_display = ref_base + indel_seq
        alt_display = ref_base
        var_type = "DEL"
    
    # Zygosity
    zygosity = determine_zygosity(af, true_gene_cn)
    
    # Quality score
    quality_score = calculate_variant_quality(indel_count, total_depth, af)
    
    # Basic filtering for indels
    filters = []
    if indel_count < MIN_ALT_READS:
        filters.append("LowAltReads")
    if af < MIN_ALT_FRACTION:
        filters.append("LowAF")
    if total_depth < MIN_TOTAL_DEPTH:
        filters.append("LowDepth")
    
    # Check homology region
    in_homology = position in homology_sites if homology_sites else False
    if in_homology:
        filters.append("HomologyRegion")
    
    filter_status = ";".join(filters) if filters else "PASS"
    
    variant_id = f"{chrom}:{position}:{ref_display}>{alt_display}"
    
    return OtoaVariantCall(
        variant_id=variant_id,
        chromosome=chrom,
        position=position,
        ref_base=ref_display,
        alt_base=alt_display,
        variant_type=var_type,
        allele_fraction=round(af, 4),
        total_depth=total_depth,
        ref_reads=total_depth - indel_count,
        alt_reads=indel_count,
        zygosity=zygosity,
        filter_status=filter_status,
        quality_score=round(quality_score, 1),
        strand_bias_pval=None,
        in_homology_region=in_homology,
        pseudogene_evidence=None,
        haplotype=None,  # Set later by phase_variants()
        phase_set=None,  # Set later by phase_variants()
    )


# =============================================================================
# Quality and Filter Functions
# =============================================================================

def calculate_strand_bias(forward_alt: int, reverse_alt: int, ref_count: int) -> Optional[float]:
    """
    Calculate strand bias using Fisher's exact test.
    """
    if forward_alt + reverse_alt < 3:
        return None
    
    # Assume reference is evenly distributed
    forward_ref = ref_count // 2
    reverse_ref = ref_count - forward_ref
    
    try:
        _, pvalue = fisher_exact([
            [forward_alt, reverse_alt],
            [max(1, forward_ref), max(1, reverse_ref)],
        ])
        return pvalue
    except Exception:
        return None


def calculate_variant_quality(alt_count: int, total_depth: int, af: float) -> float:
    """
    Calculate Phred-scaled quality score for variant.
    
    Based on binomial probability of observing alt_count reads
    if the true frequency were due to error alone.
    """
    # Assume error rate of 0.01
    error_rate = 0.01
    
    try:
        # Use scipy.stats for binomial test
        from scipy.stats import binom
        
        # Probability of seeing this many or more alt reads by chance
        pval = 1 - binom.cdf(alt_count - 1, total_depth, error_rate)
        
        # Convert to Phred scale
        if pval > 0:
            qual = -10 * math.log10(pval)
            return min(qual, 999)  # Cap at 999
        else:
            return 999
    except Exception:
        return 0


def determine_zygosity(af: float, copy_number: int) -> str:
    """
    Determine zygosity based on allele fraction and copy number.
    """
    if copy_number == 0:
        return "NO_COPY"
    elif copy_number == 1:
        # Hemizygous - single copy
        if af >= 0.8:
            return "HEMI"
        elif af >= 0.3:
            return "HEMI_MOSAIC"
        else:
            return "LOW_AF"
    elif copy_number == 2:
        # Diploid
        if af >= HOM_AF_MIN:
            return "HOM"
        elif af >= HET_AF_MIN:
            return "HET"
        else:
            return "LOW_AF"
    else:
        # High copy number - estimate copies with variant
        copies_with_var = round(af * copy_number)
        if copies_with_var == copy_number:
            return "HOM"
        elif copies_with_var == 0:
            return "LOW_AF"
        else:
            return f"HET_{copies_with_var}/{copy_number}"


def apply_variant_filters(
    af: float,
    alt_count: int,
    total_depth: int,
    strand_bias_pval: Optional[float],
    forward_alt: int,
    reverse_alt: int,
    in_homology: bool,
    pseudogene_evidence: Optional[str],
) -> str:
    """
    Apply quality filters to variant.
    
    Returns:
        Filter status string (PASS or semicolon-separated filter names)
    """
    filters = []
    
    # Depth filters
    if total_depth < 20:
        filters.append("LowDepth")
    if alt_count < 5:
        filters.append("LowAltReads")
    
    # Strand bias filter (but check strand ratio too)
    if strand_bias_pval is not None and strand_bias_pval < STRAND_BIAS_PVALUE:
        total_alt = forward_alt + reverse_alt
        if total_alt > 0:
            minor_strand_ratio = min(forward_alt, reverse_alt) / total_alt
            if minor_strand_ratio < MIN_STRAND_RATIO:
                filters.append("StrandBias")
    
    # Homology region warning
    if in_homology:
        filters.append("HomologyRegion")
    
    # Pseudogene contamination
    if pseudogene_evidence and "contamination" in pseudogene_evidence.lower():
        filters.append("PseudogeneContamination")
    
    # Borderline AF
    if HET_AF_MIN <= af < 0.25:
        filters.append("BorderlineAF")
    
    return ";".join(filters) if filters else "PASS"


def check_pseudogene_signal(
    bamfile: pysam.AlignmentFile,
    chrom: str,
    position: int,
    alt_base: str,
    pseudogene_region: Dict,
) -> Optional[str]:
    """
    Check if variant might be due to pseudogene read contamination.
    
    Strategy: If the alt_base is the predominant base at the corresponding
    pseudogene position, this suggests the variant reads may be mismapped
    from pseudogene.
    
    Args:
        bamfile: Open BAM file
        chrom: Chromosome
        position: True gene position
        alt_base: Alternate base observed
        pseudogene_region: Pseudogene region coordinates
        
    Returns:
        String describing evidence, or None
    """
    # Calculate corresponding pseudogene position
    # This requires knowing the alignment/offset between true gene and pseudogene
    # For now, we'll use a simple heuristic based on the region offset
    
    true_gene_start = pseudogene_region.get('true_gene_start', 0)
    pseudo_start = pseudogene_region.get('start', 0)
    
    if true_gene_start == 0:
        return None
    
    # Estimate corresponding position in pseudogene
    offset = position - true_gene_start
    pseudo_pos = pseudo_start + offset
    
    try:
        # Check what's at the pseudogene position
        base_counts = defaultdict(int)
        
        for pileup_column in bamfile.pileup(
            chrom,
            pseudo_pos - 1,
            pseudo_pos,
            truncate=True,
            min_mapping_quality=MIN_MAPPING_QUALITY,
        ):
            if pileup_column.pos == pseudo_pos - 1:
                for pileup_read in pileup_column.pileups:
                    if pileup_read.is_del or pileup_read.is_refskip:
                        continue
                    try:
                        base = pileup_read.alignment.query_sequence[pileup_read.query_position].upper()
                        base_counts[base] += 1
                    except (TypeError, IndexError):
                        continue
        
        total = sum(base_counts.values())
        if total > 0 and alt_base.upper() in base_counts:
            pseudo_fraction = base_counts[alt_base.upper()] / total
            if pseudo_fraction > 0.8:
                return f"Likely pseudogene contamination (pseudo has {alt_base} at {pseudo_fraction:.0%})"
        
        return None
        
    except Exception as e:
        logging.debug(f"Error checking pseudogene: {e}")
        return None


# =============================================================================
# Utility Functions
# =============================================================================

def load_homology_sites(snp_file: str) -> Set[int]:
    """
    Load positions known to be in high-homology regions from SNP file.
    
    These are positions where true gene and pseudogene sequences differ,
    which can help identify homology regions.
    
    Args:
        snp_file: Path to SNP definition file
        
    Returns:
        Set of genomic positions
    """
    homology_sites = set()
    
    try:
        with open(snp_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) >= 2:
                    try:
                        pos = int(fields[1])
                        homology_sites.add(pos)
                    except ValueError:
                        continue
    except Exception as e:
        logging.warning(f"Could not load homology sites: {e}")
    
    return homology_sites


def filter_variants_by_quality(
    variants: List[OtoaVariantCall],
    min_quality: float = 30,
    require_pass: bool = False,
) -> List[OtoaVariantCall]:
    """
    Filter variants by quality metrics.
    
    Args:
        variants: List of variant calls
        min_quality: Minimum quality score
        require_pass: If True, only return PASS variants
        
    Returns:
        Filtered list of variants
    """
    filtered = []
    
    for var in variants:
        if var.quality_score < min_quality:
            continue
        if require_pass and var.filter_status != "PASS":
            continue
        filtered.append(var)
    
    return filtered


def summarize_variants(variants: List[OtoaVariantCall]) -> Dict:
    """
    Generate summary statistics for variant calls.
    """
    if not variants:
        return {
            "total_variants": 0,
            "pass_variants": 0,
            "snvs": 0,
            "insertions": 0,
            "deletions": 0,
            "het_variants": 0,
            "hom_variants": 0,
        }
    
    pass_vars = [v for v in variants if v.filter_status == "PASS"]
    snvs = [v for v in variants if v.variant_type == "SNV"]
    insertions = [v for v in variants if v.variant_type == "INS"]
    deletions = [v for v in variants if v.variant_type == "DEL"]
    het_vars = [v for v in variants if "HET" in v.zygosity]
    hom_vars = [v for v in variants if v.zygosity == "HOM"]
    
    return {
        "total_variants": len(variants),
        "pass_variants": len(pass_vars),
        "snvs": len(snvs),
        "insertions": len(insertions),
        "deletions": len(deletions),
        "het_variants": len(het_vars),
        "hom_variants": len(hom_vars),
        "mean_af": sum(v.allele_fraction for v in variants) / len(variants),
        "mean_depth": sum(v.total_depth for v in variants) / len(variants),
        "mean_quality": sum(v.quality_score for v in variants) / len(variants),
    }


def update_variants_with_phasing(
    variants: List[OtoaVariantCall],
    phasing_result: Dict,
) -> List[OtoaVariantCall]:
    """
    Update variant calls with haplotype phasing information.

    Args:
        variants: List of OtoaVariantCall objects
        phasing_result: Output from phase_variants() containing:
            - variant_haplotypes: {variant_id: "H1"|"H2"|...}
            - phase_sets: [[var1, var2, ...], ...]

    Returns:
        List of updated OtoaVariantCall objects with haplotype info
    """
    variant_haplotypes = phasing_result.get("variant_haplotypes", {})
    phase_sets = phasing_result.get("phase_sets", [])

    # Create phase set mapping
    phase_set_map = {}
    for i, ps in enumerate(phase_sets):
        for var_id in ps:
            phase_set_map[var_id] = f"PS{i + 1}"

    updated_variants = []
    for v in variants:
        haplotype = variant_haplotypes.get(v.variant_id)
        phase_set = phase_set_map.get(v.variant_id)

        # Create new variant with updated fields
        updated = OtoaVariantCall(
            variant_id=v.variant_id,
            chromosome=v.chromosome,
            position=v.position,
            ref_base=v.ref_base,
            alt_base=v.alt_base,
            variant_type=v.variant_type,
            allele_fraction=v.allele_fraction,
            total_depth=v.total_depth,
            ref_reads=v.ref_reads,
            alt_reads=v.alt_reads,
            zygosity=v.zygosity,
            filter_status=v.filter_status,
            quality_score=v.quality_score,
            strand_bias_pval=v.strand_bias_pval,
            in_homology_region=v.in_homology_region,
            pseudogene_evidence=v.pseudogene_evidence,
            haplotype=haplotype,
            phase_set=phase_set,
        )
        updated_variants.append(updated)

    return updated_variants


def variants_to_dict(variants: List[OtoaVariantCall]) -> List[Dict]:
    """
    Convert variant calls to list of dictionaries for JSON output.
    """
    return [
        {
            "id": v.variant_id,
            "chr": v.chromosome,
            "pos": v.position,
            "ref": v.ref_base,
            "alt": v.alt_base,
            "type": v.variant_type,
            "af": v.allele_fraction,
            "depth": v.total_depth,
            "alt_reads": v.alt_reads,
            "zygosity": v.zygosity,
            "filter": v.filter_status,
            "quality": v.quality_score,
            "in_homology": v.in_homology_region,
            "pseudogene_flag": v.pseudogene_evidence,
            "haplotype": v.haplotype,
            "phase_set": v.phase_set,
        }
        for v in variants
    ]
