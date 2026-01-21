#!/usr/bin/env python3
#
# CNV and SV detection for OTOA and its pseudogene from WGS
# Author: Yu-Hui Lin <yhlin.md05@nycu.edu.tw>
# OTOAlyzer is based on CNVPanelizer (https://github.com/biostuff/CNVPanelizer) and Cyrius (https://github.com/illumina/Cyrius)
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


# --- CNVPanelizer integration module ---
# Replaces Cyrius GMM-based CNV calling with CNVPanelizer's reference-panel normalization approach.
# Executes external R script and parses output to determine total_cn and otoa for WGS data. 


import os
import subprocess
import pandas as pd
import math
import logging


def get_cn_from_cnvpanelizer(bam_file, r_script_path, output_dir, bed_file_path, reference_dir_path, reference_fasta=None):
    """
    Runs the CNVPanelizer R script and calculates total_cn and otoa (true gene CN).

    Args:
        bam_file (str): Path to the sample BAM file.
        r_script_path (str): Path to the run_CNVPanelizer.R script.
        output_dir (str): Directory where the CNVPanelizer report will be saved.
        bed_file_path (str): Path to the BED file for CNVPanelizer.
        reference_dir_path (str): Path to the directory of reference BAMs.
        reference_fasta (str, optional): Path to reference FASTA for CRAM support.

    Returns:
        tuple: A tuple containing (total_cn, otoa), or (None, None) on failure.
               - total_cn: Total copies of OTOA + OTOAP1 (from OTOA_full + OTOAP1_full)
               - otoa: Copy number of true OTOA gene only (from OTOA_unique)

               Pseudogene CN can be derived as: pseudogene_cn = total_cn - otoa
    """
    # --- Validate input files and create output directory. ---
    if not os.path.exists(bam_file):
        logging.error(f"BAM file not found: {bam_file}")
        return None, None

    if not os.path.exists(r_script_path):
        logging.error(f"R script not found: {r_script_path}")
        return None, None

    os.makedirs(output_dir, exist_ok=True)

    bam_basename = os.path.splitext(os.path.basename(bam_file))[0]
    report_filename = os.path.join(output_dir, f"{bam_basename}_CNV_exon_level_report.csv")

    # --- STEP 1: Execute CNVPanelizer R script via subprocess ---
    # Calls external R script with sample BAM, BED file, reference BAMs directory, and output directory.
    # CNVPanelizer normalizes sample depth against reference panel to calculate copy number ratios.
    logging.info(f"Running CNVPanelizer for {bam_basename}...")
    try:
        r_command = [
            "Rscript",
            r_script_path,
            bam_file,
            bed_file_path,
            reference_dir_path,
            output_dir,
        ]
        # Add reference FASTA for CRAM support if provided
        if reference_fasta:
            r_command.append(reference_fasta)
        result = subprocess.run(
            r_command, 
            check=True, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, 
            universal_newlines=True
        )
        logging.info("CNVPanelizer R script executed successfully.")
        logging.debug(f"R stdout: {result.stdout}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error executing CNVPanelizer R script for {bam_basename}:")
        logging.error(f"stderr: {e.stderr}")
        return None, None
    except FileNotFoundError:
        logging.error("Rscript not found. Please ensure R is installed and in PATH.")
        return None, None

    # --- STEP 2: Load CNVPanelizer CSV output report ---
    # Report contains MeanRatio values for each genomic region (gene-level CN estimates). 
    if not os.path.exists(report_filename):
        logging.error(f"CNVPanelizer report not found: {report_filename}")
        return None, None

    try:
        # First column (gene names) used as index
        report_df = pd.read_csv(report_filename, index_col=0)
        logging.debug(f"CNVPanelizer report loaded with {len(report_df)} regions")
    except Exception as e:
        logging.error(f"Failed to read or parse the CNVPanelizer report: {e}")
        return None, None

    # --- STEP 3: Extract MeanRatio values for key genomic regions. ---
    # Three regions extracted: OTOA_full, OTOAP1_full, and OTOA_unique (each has expected haploid CN). 
    try:
        mean_ratio_true = report_df.loc["OTOA_full", "MeanRatio"]
        mean_ratio_pseudo = report_df.loc["OTOAP1_full", "MeanRatio"]
        mean_ratio_otoa = report_df.loc["OTOA_unique", "MeanRatio"]
        logging.info(f"MeanRatio values: OTOA_full={mean_ratio_true:.3f}, "
                     f"OTOAP1_full={mean_ratio_pseudo:.3f}, OTOA_unique={mean_ratio_otoa:.3f}")
    except KeyError as e:
        logging.error(f"Could not find required gene name in the report's index: {e}")
        logging.error(f"Available regions: {list(report_df.index)}")
        return None, None

    # --- STEP 4: Calculate total_cn (OTOA_full+OTOAP1_full) and OTOA (OTOA_unique) from MeanRatio values. ---
    # Formula: total_cn = round(true_ratio*2) + round(pseudo_ratio*2)
    # Combines OTOA and OTOAP1 CN estimates to get total gene cluster copy number.
    # otoa = round(OTOA_unique*2) gives OTOA copy number.
    # Epsilon values (currently 0) allow tuning/hysteresis based on reference panel characteristics. 
    true_eps = 0  # Hysteresis window: tune on reference samples if needed
    pseudo_eps = 0  # Hysteresis window: tune on reference samples if needed
    total_cn = int(round((mean_ratio_true + true_eps) * 2) + round((mean_ratio_pseudo + pseudo_eps) * 2))

    otoa_eps = 0  # Hysteresis window: tune on reference samples if needed
    otoa = int(round((mean_ratio_otoa + otoa_eps) * 2))

    logging.info(f"Calculated total_cn: {total_cn}")
    logging.info(f"Calculated otoa (true gene CN): {otoa}")

    return total_cn, otoa
