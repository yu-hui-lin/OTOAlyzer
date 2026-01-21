#!/usr/bin/env python3
#
# CNV calling for OTOA using cnvkit flat reference approach
# Author: Yu-Hui Lin <yhlin.md05@nycu.edu.tw>
# OTOAlyzer - cnvkit integration module
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

"""
CNVkit flat reference integration for OTOA copy number calling.

This module replaces CNVPanelizer with cnvkit's flat reference approach:
- No need for a panel of normal samples
- Uses genome-wide bins (antitargets) for depth normalization
- Applies GC content correction from reference FASTA

Key advantage: Single-sample analysis without requiring matched normals or a reference panel.

cnvkit workflow:
1. Prepare targets (OTOA regions from BED file)
2. Generate antitargets (genome-wide bins for normalization)
3. Calculate coverage for sample BAM (both targets and antitargets)
4. Create flat reference with GC correction
5. Fix/normalize sample coverage against flat reference
6. Call integer copy numbers

Reference: https://cnvkit.readthedocs.io/en/stable/
"""

import os
import subprocess
import tempfile
import logging
import shutil
from typing import Tuple, Optional, Dict, List
import pandas as pd


# =============================================================================
# Constants
# =============================================================================

# cnvkit default parameters
DEFAULT_TARGET_AVG_SIZE = 50000   # Average target bin size (bp)
DEFAULT_ANTITARGET_AVG_SIZE = 100000  # Average antitarget bin size (bp)
DEFAULT_ANTITARGET_MIN_SIZE = 10000   # Minimum antitarget bin size (bp)

# Copy number conversion
DIPLOID_LOG2 = 0.0  # log2 ratio for diploid (CN=2)


# =============================================================================
# Utility Functions
# =============================================================================

def check_cnvkit_installed() -> bool:
    """
    Check if cnvkit is installed and accessible.

    Returns:
        True if cnvkit is available, False otherwise
    """
    try:
        result = subprocess.run(
            ["cnvkit.py", "version"],
            capture_output=True,
            text=True,
            check=False
        )
        if result.returncode == 0:
            version = result.stdout.strip() or result.stderr.strip()
            logging.debug(f"cnvkit version: {version}")
            return True
        return False
    except FileNotFoundError:
        return False


def run_cnvkit_command(args: List[str], description: str) -> subprocess.CompletedProcess:
    """
    Run a cnvkit command with error handling.

    Args:
        args: Command arguments (without 'cnvkit.py' prefix)
        description: Description of the command for logging

    Returns:
        CompletedProcess object

    Raises:
        RuntimeError: If command fails
    """
    cmd = ["cnvkit.py"] + args
    logging.debug(f"Running: {' '.join(cmd)}")

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        logging.debug(f"{description} completed successfully")
        return result
    except subprocess.CalledProcessError as e:
        logging.error(f"{description} failed:")
        logging.error(f"  Command: {' '.join(cmd)}")
        logging.error(f"  stdout: {e.stdout}")
        logging.error(f"  stderr: {e.stderr}")
        raise RuntimeError(f"{description} failed: {e.stderr}") from e


def log2_to_copy_number(log2_ratio: float, ploidy: int = 2) -> int:
    """
    Convert log2 ratio to integer copy number.

    Formula: CN = round(ploidy * 2^log2_ratio)

    Args:
        log2_ratio: Log2 ratio from cnvkit
        ploidy: Expected ploidy (default 2 for diploid)

    Returns:
        Integer copy number
    """
    cn = ploidy * (2 ** log2_ratio)
    return max(0, round(cn))


# =============================================================================
# CNVkit Workflow Functions
# =============================================================================

def prepare_targets(
    bed_file: str,
    output_dir: str,
    reference_fasta: str,
    avg_size: int = DEFAULT_TARGET_AVG_SIZE
) -> str:
    """
    Prepare target regions from BED file with GC annotation.

    Args:
        bed_file: Input BED file with OTOA regions
        output_dir: Directory for output files
        reference_fasta: Path to reference FASTA
        avg_size: Average target bin size

    Returns:
        Path to prepared target file (.target.bed)
    """
    output_file = os.path.join(output_dir, "otoa.target.bed")

    run_cnvkit_command(
        ["target", bed_file,
         #"--fasta", reference_fasta,
         "--avg-size", str(avg_size),
         "-o", output_file],
        "Target preparation"
    )

    return output_file


def prepare_antitargets(
    target_file: str,
    output_dir: str,
    reference_fasta: str,
    access_file: Optional[str] = None,
    avg_size: int = DEFAULT_ANTITARGET_AVG_SIZE,
    min_size: int = DEFAULT_ANTITARGET_MIN_SIZE
) -> str:
    """
    Generate antitarget regions for genome-wide normalization.

    Antitargets are regions outside the target BED that provide
    background coverage for normalization.

    Args:
        target_file: Prepared target BED file
        output_dir: Directory for output files
        reference_fasta: Path to reference FASTA
        access_file: Optional accessible regions file (to exclude problematic regions)
        avg_size: Average antitarget bin size
        min_size: Minimum antitarget bin size

    Returns:
        Path to antitarget file (.antitarget.bed)
    """
    output_file = os.path.join(output_dir, "otoa.antitarget.bed")

    args = [
        "antitarget", target_file,
        "--avg-size", str(avg_size),
        "--min-size", str(min_size),
        "-o", output_file
    ]

    if access_file:
        args.extend(["--access", access_file])

    run_cnvkit_command(args, "Antitarget generation")

    return output_file


def calculate_coverage(
    bam_file: str,
    bed_file: str,
    output_file: str,
    reference_fasta: Optional[str] = None
) -> str:
    """
    Calculate coverage for a BAM file over specified regions.

    Args:
        bam_file: Input BAM/CRAM file
        bed_file: Target or antitarget BED file
        output_file: Output coverage file (.targetcoverage.cnn or .antitargetcoverage.cnn)
        reference_fasta: Reference FASTA (required for CRAM files)

    Returns:
        Path to coverage file
    """
    args = ["coverage", bam_file, bed_file, "-o", output_file]

    if reference_fasta:
        args.extend(["--fasta", reference_fasta])

    run_cnvkit_command(args, f"Coverage calculation for {os.path.basename(bed_file)}")

    return output_file


def create_flat_reference(
    target_file: str,
    antitarget_file: str,
    output_file: str,
    reference_fasta: str
) -> str:
    """
    Create a flat reference for single-sample analysis.

    A flat reference assumes neutral copy number (log2 = 0) for all regions,
    with GC content calculated from the reference FASTA for correction.

    Args:
        target_file: Target BED file
        antitarget_file: Antitarget BED file
        output_file: Output reference file (.cnn)
        reference_fasta: Path to reference FASTA

    Returns:
        Path to flat reference file
    """
    run_cnvkit_command(
        ["reference",
         "--fasta", reference_fasta,
         "--targets", target_file,
         "--antitargets", antitarget_file,
         "-o", output_file],
        "Flat reference creation"
    )

    return output_file


def fix_coverage(
    target_coverage: str,
    antitarget_coverage: str,
    reference_file: str,
    output_file: str
) -> str:
    """
    Normalize sample coverage against the flat reference.

    This step applies GC correction and normalizes using antitarget coverage.

    Args:
        target_coverage: Sample target coverage file (.targetcoverage.cnn)
        antitarget_coverage: Sample antitarget coverage file (.antitargetcoverage.cnn)
        reference_file: Flat reference file (.cnn)
        output_file: Output normalized ratios file (.cnr)

    Returns:
        Path to fixed/normalized ratios file
    """
    run_cnvkit_command(
        ["fix", target_coverage, antitarget_coverage, reference_file,
         "-o", output_file],
        "Coverage normalization"
    )

    return output_file


def call_copy_numbers(
    cnr_file: str,
    output_file: str,
    ploidy: int = 2,
    method: str = "threshold"
) -> str:
    """
    Call integer copy numbers from normalized ratios.

    Args:
        cnr_file: Normalized ratios file (.cnr)
        output_file: Output calls file (.cns or .call.cns)
        ploidy: Expected sample ploidy (default 2)
        method: Calling method ('threshold', 'clonal', 'none')

    Returns:
        Path to calls file
    """
    run_cnvkit_command(
        ["call", cnr_file,
         "--ploidy", str(ploidy),
         "--method", method,
         "-o", output_file],
        "Copy number calling"
    )

    return output_file


def segment_ratios(
    cnr_file: str,
    output_file: str,
    method: str = "cbs"
) -> str:
    """
    Segment normalized ratios to smooth noisy data.

    Args:
        cnr_file: Normalized ratios file (.cnr)
        output_file: Output segments file (.cns)
        method: Segmentation method ('cbs', 'flasso', 'haar', 'none')

    Returns:
        Path to segments file
    """
    run_cnvkit_command(
        ["segment", cnr_file,
         "--method", method,
         "-o", output_file],
        "Segmentation"
    )

    return output_file


# =============================================================================
# Main Interface Function
# =============================================================================

def get_cn_from_cnvkit(
    bam_file: str,
    bed_file_path: str,
    reference_fasta: str,
    output_dir: str,
    access_file: Optional[str] = None,
    keep_intermediate: bool = False
) -> Tuple[Optional[int], Optional[int]]:
    """
    Determine OTOA copy numbers using cnvkit flat reference approach.

    This function replaces get_cn_from_cnvpanelizer() with a single-sample
    analysis that uses genome-wide bins for normalization instead of a
    panel of normal samples.

    Args:
        bam_file: Path to sample BAM/CRAM file
        bed_file_path: Path to BED file with OTOA regions
            Expected regions: OTOA_full, OTOAP1_full, OTOA_unique
        reference_fasta: Path to reference FASTA (REQUIRED for GC correction)
        output_dir: Directory for output files
        access_file: Optional accessible regions file
        keep_intermediate: If True, keep all intermediate files

    Returns:
        Tuple of (total_cn, otoa_cn):
        - total_cn: Total copies of OTOA + OTOAP1 (from OTOA_full + OTOAP1_full)
        - otoa_cn: Copy number of true OTOA gene only (from OTOA_unique)

        Returns (None, None) on failure.

        Pseudogene CN can be derived as: pseudogene_cn = total_cn - otoa_cn
    """
    # --- Validate inputs ---
    if not os.path.exists(bam_file):
        logging.error(f"BAM file not found: {bam_file}")
        return None, None

    if not os.path.exists(bed_file_path):
        logging.error(f"BED file not found: {bed_file_path}")
        return None, None

    if not os.path.exists(reference_fasta):
        logging.error(f"Reference FASTA not found: {reference_fasta}")
        logging.error("Reference FASTA is REQUIRED for cnvkit GC correction")
        return None, None

    if not check_cnvkit_installed():
        logging.error("cnvkit is not installed or not in PATH")
        logging.error("Install with: pip install cnvkit")
        return None, None

    bam_basename = os.path.splitext(os.path.basename(bam_file))[0]
    logging.info(f"Running cnvkit flat reference analysis for {bam_basename}...")

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Use temporary directory for intermediate files if not keeping them
    if keep_intermediate:
        work_dir = output_dir
    else:
        work_dir = tempfile.mkdtemp(prefix="cnvkit_otoa_")

    try:
        # =====================================================================
        # STEP 1: Prepare target and antitarget regions
        # =====================================================================
        logging.info("STEP 1: Preparing target and antitarget regions...")

        target_file = prepare_targets(
            bed_file=bed_file_path,
            output_dir=work_dir,
            reference_fasta=reference_fasta
        )

        antitarget_file = prepare_antitargets(
            target_file=target_file,
            output_dir=work_dir,
            reference_fasta=reference_fasta,
            access_file=access_file
        )

        # =====================================================================
        # STEP 2: Calculate coverage
        # =====================================================================
        logging.info("STEP 2: Calculating coverage for sample...")

        target_cov_file = os.path.join(work_dir, f"{bam_basename}.targetcoverage.cnn")
        antitarget_cov_file = os.path.join(work_dir, f"{bam_basename}.antitargetcoverage.cnn")

        calculate_coverage(
            bam_file=bam_file,
            bed_file=target_file,
            output_file=target_cov_file,
            reference_fasta=reference_fasta
        )

        calculate_coverage(
            bam_file=bam_file,
            bed_file=antitarget_file,
            output_file=antitarget_cov_file,
            reference_fasta=reference_fasta
        )

        # =====================================================================
        # STEP 3: Create flat reference
        # =====================================================================
        logging.info("STEP 3: Creating flat reference with GC correction...")

        flat_ref_file = os.path.join(work_dir, "flat_reference.cnn")

        create_flat_reference(
            target_file=target_file,
            antitarget_file=antitarget_file,
            output_file=flat_ref_file,
            reference_fasta=reference_fasta
        )

        # =====================================================================
        # STEP 4: Fix/Normalize coverage
        # =====================================================================
        logging.info("STEP 4: Normalizing coverage against flat reference...")

        cnr_file = os.path.join(work_dir, f"{bam_basename}.cnr")

        fix_coverage(
            target_coverage=target_cov_file,
            antitarget_coverage=antitarget_cov_file,
            reference_file=flat_ref_file,
            output_file=cnr_file
        )

        # =====================================================================
        # STEP 5: Call copy numbers
        # =====================================================================
        logging.info("STEP 5: Calling copy numbers...")

        call_file = os.path.join(output_dir, f"{bam_basename}.call.cns")

        call_copy_numbers(
            cnr_file=cnr_file,
            output_file=call_file
        )

        # =====================================================================
        # STEP 6: Parse results and extract OTOA CN values
        # =====================================================================
        logging.info("STEP 6: Parsing cnvkit results...")

        total_cn, otoa_cn = parse_cnvkit_results(cnr_file, call_file)

        if total_cn is not None:
            logging.info(f"cnvkit results:")
            logging.info(f"  Total CN (OTOA+OTOAP1): {total_cn}")
            logging.info(f"  True gene CN (OTOA_unique): {otoa_cn}")

        return total_cn, otoa_cn

    except Exception as e:
        logging.error(f"cnvkit analysis failed: {e}")
        return None, None

    finally:
        # Clean up temporary directory if not keeping intermediate files
        if not keep_intermediate and work_dir != output_dir:
            try:
                shutil.rmtree(work_dir)
            except Exception as e:
                logging.warning(f"Failed to clean up temporary directory: {e}")


def parse_cnvkit_results(
    cnr_file: str,
    call_file: str
) -> Tuple[Optional[int], Optional[int]]:
    """
    Parse cnvkit output files to extract OTOA copy numbers.

    Looks for regions named:
    - OTOA_full: Full OTOA region
    - OTOAP1_full: Full pseudogene region
    - OTOA_unique: OTOA-unique region (for true gene CN)

    Args:
        cnr_file: Path to .cnr file (normalized ratios)
        call_file: Path to .call.cns file (called copy numbers)

    Returns:
        Tuple of (total_cn, otoa_cn)
    """
    try:
        # Read the CNR file for log2 ratios
        cnr_df = pd.read_csv(cnr_file, sep='\t')
        logging.debug(f"CNR file loaded with {len(cnr_df)} bins")
        logging.debug(f"CNR columns: {list(cnr_df.columns)}")

        # Read the call file for integer copy numbers
        call_df = pd.read_csv(call_file, sep='\t')
        logging.debug(f"Call file loaded with {len(call_df)} segments")

        # Look for OTOA regions in the gene column
        # cnvkit may use 'gene' or 'name' column depending on version
        gene_col = 'gene' if 'gene' in cnr_df.columns else 'name' if 'name' in cnr_df.columns else None

        if gene_col is None:
            logging.error(f"Could not find gene/name column in CNR file. Columns: {list(cnr_df.columns)}")
            # Try to match by coordinates
            return parse_by_coordinates(cnr_df, call_df)

        # Extract log2 ratios for each region
        otoa_full_log2 = None
        otoap1_full_log2 = None
        otoa_unique_log2 = None

        for idx, row in cnr_df.iterrows():
            gene_name = str(row[gene_col]) if pd.notna(row[gene_col]) else ""
            if 'OTOA_full' in gene_name:
                otoa_full_log2 = row['log2']
                logging.debug(f"OTOA_full log2: {otoa_full_log2}")
            elif 'OTOAP1_full' in gene_name:
                otoap1_full_log2 = row['log2']
                logging.debug(f"OTOAP1_full log2: {otoap1_full_log2}")
            elif 'OTOA_unique' in gene_name:
                otoa_unique_log2 = row['log2']
                logging.debug(f"OTOA_unique log2: {otoa_unique_log2}")

        # Calculate copy numbers from log2 ratios
        if otoa_full_log2 is None or otoap1_full_log2 is None:
            logging.error("Could not find OTOA_full or OTOAP1_full regions in CNR file")
            logging.debug(f"Available genes: {cnr_df[gene_col].unique()[:20]}")
            return None, None

        # Convert log2 to CN
        # Formula: CN = 2 * 2^log2 (for diploid baseline)
        cn_otoa_full = log2_to_copy_number(otoa_full_log2)
        cn_otoap1_full = log2_to_copy_number(otoap1_full_log2)

        # Total CN = OTOA + OTOAP1
        total_cn = cn_otoa_full + cn_otoap1_full

        # True gene CN from OTOA_unique if available, otherwise use OTOA_full
        if otoa_unique_log2 is not None:
            otoa_cn = log2_to_copy_number(otoa_unique_log2)
        else:
            logging.warning("OTOA_unique region not found, using OTOA_full for true gene CN")
            otoa_cn = cn_otoa_full

        logging.info(f"Log2 ratios: OTOA_full={otoa_full_log2}, "
                     f"OTOAP1_full={otoap1_full_log2}, "
                     f"OTOA_unique={otoa_unique_log2 if otoa_unique_log2 else 'N/A'}")

        return total_cn, otoa_cn

    except Exception as e:
        logging.error(f"Failed to parse cnvkit results: {e}")
        import traceback
        traceback.print_exc()
        return None, None


def parse_by_coordinates(
    cnr_df: pd.DataFrame,
    call_df: pd.DataFrame
) -> Tuple[Optional[int], Optional[int]]:
    """
    Parse cnvkit results by matching coordinates instead of gene names.

    Expected OTOA regions (GRCh38):
    - OTOA_full: chr16:21,677,882-21,761,493
    - OTOAP1_full: chr16:21,493,807-21,577,518
    - OTOA_unique: chr16:21,729,000-21,760,000

    Args:
        cnr_df: DataFrame from .cnr file
        call_df: DataFrame from .call.cns file

    Returns:
        Tuple of (total_cn, otoa_cn)
    """
    # Define expected regions (GRCh38 coordinates)
    expected_regions = {
        'OTOA_full': ('chr16', 21677882, 21761493),
        'OTOAP1_full': ('chr16', 21493807, 21577518),
        'OTOA_unique': ('chr16', 21729000, 21760000),
    }

    log2_values = {}

    for region_name, (chrom, start, end) in expected_regions.items():
        # Find matching bin by overlapping coordinates
        mask = (
            (cnr_df['chromosome'] == chrom) &
            (cnr_df['start'] >= start - 1000) &
            (cnr_df['end'] <= end + 1000)
        )

        matches = cnr_df[mask]
        if len(matches) > 0:
            # Use weighted average if multiple bins
            log2_values[region_name] = matches['log2'].mean()
            logging.debug(f"{region_name} log2 (by coords): {log2_values[region_name]:.3f}")
        else:
            logging.warning(f"Could not find {region_name} region by coordinates")

    if 'OTOA_full' not in log2_values or 'OTOAP1_full' not in log2_values:
        logging.error("Could not find required OTOA regions by coordinates")
        return None, None

    cn_otoa_full = log2_to_copy_number(log2_values['OTOA_full'])
    cn_otoap1_full = log2_to_copy_number(log2_values['OTOAP1_full'])

    total_cn = cn_otoa_full + cn_otoap1_full

    if 'OTOA_unique' in log2_values:
        otoa_cn = log2_to_copy_number(log2_values['OTOA_unique'])
    else:
        otoa_cn = cn_otoa_full

    return total_cn, otoa_cn


# =============================================================================
# Alternative: Direct CNR Parsing Without Full Pipeline
# =============================================================================

def get_cn_from_cnvkit_simple(
    bam_file: str,
    bed_file_path: str,
    reference_fasta: str,
    output_dir: str
) -> Tuple[Optional[int], Optional[int]]:
    """
    Simplified cnvkit analysis using batch mode.

    This is an alternative approach that uses cnvkit's batch command
    with the --normal flag (no normal samples) to create a flat reference
    and process the sample in one command.

    Args:
        bam_file: Path to sample BAM/CRAM file
        bed_file_path: Path to BED file with OTOA regions
        reference_fasta: Path to reference FASTA
        output_dir: Directory for output files

    Returns:
        Tuple of (total_cn, otoa_cn)
    """
    if not check_cnvkit_installed():
        logging.error("cnvkit is not installed or not in PATH")
        return None, None

    bam_basename = os.path.splitext(os.path.basename(bam_file))[0]

    os.makedirs(output_dir, exist_ok=True)

    try:
        # Use batch command with --normal (no normal samples = flat reference)
        run_cnvkit_command(
            ["batch", bam_file,
             "--normal",  # Empty = flat reference
             "--targets", bed_file_path,
             "--fasta", reference_fasta,
             "--output-dir", output_dir,
             "--output-reference", os.path.join(output_dir, "flat_reference.cnn")],
            "cnvkit batch (flat reference mode)"
        )

        # Parse the output CNR file
        cnr_file = os.path.join(output_dir, f"{bam_basename}.cnr")
        call_file = os.path.join(output_dir, f"{bam_basename}.call.cns")

        if os.path.exists(cnr_file):
            return parse_cnvkit_results(cnr_file, call_file)
        else:
            logging.error(f"CNR file not found: {cnr_file}")
            return None, None

    except Exception as e:
        logging.error(f"cnvkit batch failed: {e}")
        return None, None


# =============================================================================
# Testing/Debug Functions
# =============================================================================

if __name__ == "__main__":
    import sys

    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s - %(levelname)s - %(message)s"
    )

    if len(sys.argv) < 4:
        print("Usage: python cnvkit_cn.py <bam_file> <bed_file> <reference_fasta> [output_dir]")
        print("\nThis module provides cnvkit flat reference analysis for OTOA CNV calling.")
        print("\nExample:")
        print("  python cnvkit_cn.py sample.bam OTOA_region_38.bed hg38.fa ./output")
        sys.exit(1)

    bam_file = sys.argv[1]
    bed_file = sys.argv[2]
    ref_fasta = sys.argv[3]
    output_dir = sys.argv[4] if len(sys.argv) > 4 else "./cnvkit_output"

    print(f"Testing cnvkit flat reference analysis...")
    print(f"  BAM: {bam_file}")
    print(f"  BED: {bed_file}")
    print(f"  Reference: {ref_fasta}")
    print(f"  Output: {output_dir}")

    total_cn, otoa_cn = get_cn_from_cnvkit(
        bam_file=bam_file,
        bed_file_path=bed_file,
        reference_fasta=ref_fasta,
        output_dir=output_dir,
        keep_intermediate=True
    )

    if total_cn is not None:
        print(f"\nResults:")
        print(f"  Total CN (OTOA+OTOAP1): {total_cn}")
        print(f"  True gene CN (OTOA): {otoa_cn}")
        print(f"  Pseudogene CN: {total_cn - otoa_cn}")
    else:
        print("\nAnalysis failed. Check logs for details.")
        sys.exit(1)
