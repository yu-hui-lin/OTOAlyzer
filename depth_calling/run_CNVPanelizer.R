#!/usr/bin/env Rscript
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

# --- CNVPanelizer R wrapper script ---
# This script implements reference-panel based CNV detection for OTOA and its pseudogene.
# It replaces Cyrius's GMM approach with bootstrap-based statistical analysis of normalized read counts.

# --- Personal library setup: Install packages to user's home directory to avoid permission issues on shared systems. ---
personal_lib_path <- file.path(Sys.getenv("HOME"), "R", "library")
if (!dir.exists(personal_lib_path)) {
  dir.create(personal_lib_path, recursive = TRUE)
}
.libPaths(c(personal_lib_path, .libPaths()))

# --- Parse command-line arguments: sample BAM, BED file, reference BAM directory, output directory, [reference FASTA]. ---
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4 || length(args) > 5) {
  stop("Usage: Rscript run_CNVPanelizer.R <sample_bam> <bed_file> <reference_dir> <output_dir> [reference_fasta]")
}

sampleFilenames <- args[1]
bedFilepath <- args[2]
referenceDirectory <- args[3]
output_dir <- args[4]
referenceFasta <- if (length(args) >= 5) args[5] else NULL
BAM_BASENAME <- tools::file_path_sans_ext(basename(sampleFilenames))

# --- Configure output paths. ---
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
outputReportFilename <- file.path(output_dir, paste0(BAM_BASENAME, "_CNV_exon_level_report.csv"))

# --- Main CNVPanelizer workflow with error handling. ---
tryCatch({
  
  # --- STEP 1: Install and load CNVPanelizer from Bioconductor ---
  # BiocManager handles package dependencies automatically. 
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos="https://cloud.r-project.org")
  }
  if (!requireNamespace("CNVPanelizer", quietly = TRUE)) {
    BiocManager::install("CNVPanelizer", update=FALSE)
  }
  suppressPackageStartupMessages(library(CNVPanelizer))
  message("CNVPanelizer package loaded.")
  
  # --- STEP 2: Parse BED file to extract genomic ranges and gene names ---
  # ampliconColumnNumber=4: 4th column contains gene/amplicon names.
  # split parameter set to "THIS_SHOULD_NOT_SPLIT_ANYTHING" (a nonsense string) prevents unwanted name splitting.
  # Fallback logic ensures unique identifiers for each amplicon region.
  ampliconColumnNumber <- 4
  
  message(paste("Reading BED file for CNVPanelizer:", bedFilepath))
  message(paste("Using ampliconColumnNumber:", ampliconColumnNumber))
  
  genomicRangesFromBed <- BedToGenomicRanges(
    bedFilepath,
    ampliconColumn = ampliconColumnNumber,
    split = "THIS_SHOULD_NOT_SPLIT_ANYTHING",  # A nonsense string to prevent unwanted name splitting.
    doReduce = FALSE
  )
  
  metadataFromGenomicRanges <- elementMetadata(genomicRangesFromBed)
  geneNamesPerAmplicon <- metadataFromGenomicRanges["geneNames"][, 1]
  uniqueAmpliconNames_for_functions <- metadataFromGenomicRanges["ampliconNames"][,1]
  
  # --- Fallback: Use geneNames if ampliconNames metadata not suitable. ---
  if (all(is.na(uniqueAmpliconNames_for_functions)) || length(uniqueAmpliconNames_for_functions) != length(geneNamesPerAmplicon)) {
    message("Primary 'ampliconNames' metadata was not suitable, using 'geneNamesPerAmplicon' as the unique identifiers for functions.")
    uniqueAmpliconNames_for_functions <- geneNamesPerAmplicon
  } else {
    message("Using 'ampliconNames' metadata as the unique identifiers for functions.")
  }
  
  # --- Validate uniqueness: Duplicates will cause incorrect CNV calculations. ---
  if (any(duplicated(uniqueAmpliconNames_for_functions))) {
    warning("CNVPANELIZER WARNING: The 'uniqueAmpliconNames_for_functions' derived by CNVPanelizer contain duplicates. This may lead to incorrect results.")
  } else {
    message("CNVPANELIZER INFO: 'uniqueAmpliconNames_for_functions' are unique.")
  }
  
  message(paste("Extracted", length(uniqueAmpliconNames_for_functions), "unique identifiers to be used for amplicons by CNVPanelizer."))
  
  # --- STEP 3: Load reference and sample alignment files (BAM or CRAM). ---
  # Reference panel should contain enough diploid samples for robust normalization (20+ samples recommended).
  # Fewer references may reduce CNV detection accuracy.
  # Note: CRAM files require reference FASTA set via Rsamtools (e.g., Rsamtools::setFaFile())
  referenceFilenames_bam <- list.files(path = referenceDirectory, pattern = "\\.bam$", full.names = TRUE)
  referenceFilenames_cram <- list.files(path = referenceDirectory, pattern = "\\.cram$", full.names = TRUE)
  referenceFilenames <- c(referenceFilenames_bam, referenceFilenames_cram)
  expectedReferenceCount <- 20
  if (length(referenceFilenames) != expectedReferenceCount) {
    warning(paste("Expected", expectedReferenceCount, "reference alignment files, but found", length(referenceFilenames), "in", referenceDirectory))
  } else {
    message(paste("Found", length(referenceFilenames), "reference alignment files (BAM/CRAM)."))
  }
  
  if (!file.exists(sampleFilenames)) {
    stop(paste("Test sample alignment file not found at:", sampleFilenames))
  } else {
    message(paste("Test sample alignment file:", sampleFilenames))
  }
  
  # --- Check if using CRAM files and configure reference FASTA ---
  has_cram <- grepl("\\.cram$", sampleFilenames) || any(grepl("\\.cram$", referenceFilenames))
  if (has_cram) {
    message("NOTE: CRAM files detected. CNVPanelizer uses Rsamtools internally.")
    if (!is.null(referenceFasta) && file.exists(referenceFasta)) {
      message(paste("Setting reference FASTA for CRAM support:", referenceFasta))
      suppressPackageStartupMessages(library(Rsamtools))
      Rsamtools::setFaFile(referenceFasta)
    } else if (!is.null(referenceFasta)) {
      warning(paste("Reference FASTA not found:", referenceFasta))
      message("If CRAM reading fails, ensure reference FASTA is configured:")
      message("  - Set CRAM_REFERENCE environment variable, OR")
      message("  - Provide valid path via command line argument")
    } else {
      message("No reference FASTA provided. If CRAM reading fails, ensure reference FASTA is configured:")
      message("  - Set CRAM_REFERENCE environment variable, OR")
      message("  - Pass reference FASTA path as 5th argument")
    }
  }
  
  # --- STEP 4: Count reads for each amplicon region in reference and sample BAMs. ---
  # minimumMappingQuality=20: Filter low-quality alignments.
  # removePcrDuplicates=FALSE: Keep duplicates (panel data typically has high duplication).
  minimumMappingQualityVal <- 20 
  removePcrDuplicates <- FALSE
  
  message("Counting reads for reference samples...")
  referenceReadCounts <- ReadCountsFromBam(
    bamFilenames = referenceFilenames,
    gr = genomicRangesFromBed,
    sampleNames = basename(referenceFilenames),
    ampliconNames = uniqueAmpliconNames_for_functions,
    minimumMappingQuality = minimumMappingQualityVal,
    removeDup = removePcrDuplicates
  )
  message("Finished counting reads for reference samples.")
  
  message("Counting reads for the test sample...")
  sampleReadCounts <- ReadCountsFromBam(
    bamFilenames = sampleFilenames,
    gr = genomicRangesFromBed,
    sampleNames = basename(sampleFilenames),
    ampliconNames = uniqueAmpliconNames_for_functions,
    minimumMappingQuality = minimumMappingQualityVal,
    removeDup = removePcrDuplicates
  )
  message("Finished counting reads for the test sample.")
  
  # --- Validate read count dimensions match amplicon count. ---
  if (nrow(referenceReadCounts) != length(geneNamesPerAmplicon) || nrow(sampleReadCounts) != length(geneNamesPerAmplicon)) {
    stop("Mismatch between number of amplicons and rows in read count matrices.")
  }
  
  # --- STEP 5: Normalize read counts using TMM (Trimmed Mean of M-values) method. ---
  # TMM adjusts for sequencing depth and compositional biases between samples.
  # Output: sample and reference matrices with normalized counts for each amplicon.
  message("Normalizing read counts...")
  normalizationMethodVal <- "tmm"
  normalizedReadCounts <- CombinedNormalizedCounts(
    sampleCounts = sampleReadCounts,
    referenceCounts = referenceReadCounts,
    method = normalizationMethodVal, 
    ampliconNames = uniqueAmpliconNames_for_functions
  )
  
  # --- Handle different return structures from CNVPanelizer versions. ---
  samplesNormalizedReadCounts <- normalizedReadCounts[["samples"]]
  if (is.list(samplesNormalizedReadCounts) && length(samplesNormalizedReadCounts) == 1) {
    samplesNormalizedReadCounts <- samplesNormalizedReadCounts[[1]]
  }
  
  referenceNormalizedReadCounts <- normalizedReadCounts[["reference"]]
  if (is.list(referenceNormalizedReadCounts) && length(referenceNormalizedReadCounts) == 1) {
    referenceNormalizedReadCounts <- referenceNormalizedReadCounts[[1]]
  }
  message("Normalization complete.")
  
  # --- STEP 6: Bootstrap-based CNV calling. ---
  # Generates 10,000 bootstrap replicates by resampling reference panel to estimate CN variance.
  # This provides statistical confidence for CNV calls by modeling expected variation in diploid samples.
  numReplicates <- 10000
  message(paste("Performing bootstrap-based CNV analysis with", numReplicates, "replicates..."))
  bootList <- BootList(
    geneNames = geneNamesPerAmplicon,
    sampleMatrix = samplesNormalizedReadCounts,
    refmat = referenceNormalizedReadCounts,
    replicates = numReplicates
  )
  message("Bootstrap analysis complete.")
  
  # --- STEP 7: Estimate background noise levels from reference panel. ---
  # significanceLevel=0.1: Alpha level for statistical tests.
  # robust=TRUE: Use robust statistics (median/MAD) instead of mean/SD to handle outliers.
  significanceLevelBg <- 0.1
  useRobustStats <- TRUE
  message("Estimating background noise...")
  backgroundNoise <- Background(
    geneNames = geneNamesPerAmplicon,
    samplesNormalizedReadCounts = samplesNormalizedReadCounts,
    referenceNormalizedReadCounts = referenceNormalizedReadCounts,
    bootList = bootList,
    replicates = numReplicates,
    significanceLevel = significanceLevelBg,
    robust = useRobustStats
  )
  message("Background noise estimation complete.")
  
  # --- STEP 8: Generate final CNV report with MeanRatio values. ---
  # MeanRatio: Normalized sample depth / normalized reference depth.
  # MeanRatio ≈ 1.0 indicates diploid state, 0.5 one copy deletion, 1.5 one copy duplication, etc.
  # Report includes statistical significance and confidence intervals. 
  message("Generating report tables...")
  reportTables <- ReportTables(
    geneNames = geneNamesPerAmplicon,
    samplesNormalizedReadCounts = samplesNormalizedReadCounts,
    referenceNormalizedReadCounts = referenceNormalizedReadCounts,
    bootList = bootList,
    backgroundNoise = backgroundNoise
  )
  
  cnvReportForTestSample <- reportTables[[1]]
  message("CNV Report for sample:", sampleFilenames)
  
  # --- Ensure gene names are rownames for easy lookup. ---
  if (is.null(rownames(cnvReportForTestSample))) {
    if (length(geneNamesPerAmplicon) == nrow(cnvReportForTestSample)) {
      rownames(cnvReportForTestSample) <- geneNamesPerAmplicon
    } else {
      warning("Could not set rownames for cnvReportForTestSample as length does not match.")
    }
  }
  print(cnvReportForTestSample)
  
  # --- Write CSV report for Python parser to extract MeanRatio values. ---
  write.csv(cnvReportForTestSample, file = outputReportFilename, row.names = TRUE)
  message(paste("CNV report (exon-level) saved to:", outputReportFilename))
  
}, error = function(e) {
  message("ERROR: ", e$message)
  quit(status=1)
})
