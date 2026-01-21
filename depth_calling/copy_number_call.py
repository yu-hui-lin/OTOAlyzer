#!/usr/bin/env python3
#
# CNV and SV detection for OTOA and its pseudogene from WGS
# Author: Yu-Hui Lin <yhlin.md05@nycu.edu.tw>
# OTOAlyzer is based on CNVPanelizer and Cyrius
# Original Cyrius Copyright (c) 2019-2020 Illumina, Inc.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

"""
Copy number calling module for OTOAlyzer.

Uses Poisson likelihood model to determine copy number at each SNP site
based on allele-specific read counts.
"""

from scipy.stats import poisson


# Posterior probability cutoff for confident CN calls
POSTERIOR_CUTOFF_STRINGENT = 0.9

# Error rate for Poisson model (standard WGS error rate)
ERROR_RATE = 0.01

# Threshold for resolving ambiguous CN=2 vs CN=3 calls
AMBIGUOUS_CALL_THRESHOLD = 0.85


def call_reg1_cn(full_cn, count_reg1, count_reg2, min_read=0):
    """
    Return the reg1 (OTOA true gene) copy number call at each site 
    based on Poisson likelihood, with a minimum read support cutoff.
    
    Args:
        full_cn: Total copy number (OTOA + pseudogene)
        count_reg1: Read count supporting true gene allele
        count_reg2: Read count supporting pseudogene allele
        min_read: Minimum read count required
        
    Returns:
        List with either:
        - Single value [CN] for confident calls
        - [CN1, prob1, CN2, prob2] for ambiguous calls with two candidates
        - [None] for failed calls
    """
    if full_cn is None:
        return [None]
    if full_cn == 0:
        return [0]
    if full_cn == 1 and count_reg1 > min_read:
        return [1]
    
    prob = []
    nsum = count_reg1 + count_reg2
    if nsum == 0:
        return [None]
    
    # Calculate Poisson likelihood for each possible CN
    for i in range(full_cn + 1):
        depthexpected = float(nsum) * float(i) / float(full_cn)
        if i == 0:
            depthexpected = (ERROR_RATE / 3) * float(nsum)
        if i == full_cn:
            depthexpected = float(nsum) - ERROR_RATE * float(nsum)
        if count_reg1 <= count_reg2:
            prob.append(poisson.pmf(int(count_reg1), depthexpected))
        else:
            prob.append(poisson.pmf(int(count_reg2), depthexpected))
    
    sum_prob = sum(prob)
    if sum_prob == 0:
        return [None]
    
    # Convert to posterior probabilities
    post_prob = [float(a) / float(sum_prob) for a in prob]
    if count_reg2 < count_reg1:
        post_prob = post_prob[::-1]
    
    post_prob_sorted = sorted(post_prob, reverse=True)
    
    # Handle case where CN=0 is most likely but we have reads
    if (
        post_prob.index(post_prob_sorted[0]) != 0
        and count_reg1 <= min_read
        and count_reg2 >= min_read
    ):
        return [0]
    
    # Return confident call if probability exceeds threshold
    if post_prob_sorted[0] >= POSTERIOR_CUTOFF_STRINGENT:
        return [post_prob.index(post_prob_sorted[0])]
    
    # Return two most likely scenarios for ambiguous calls
    cn_prob_filtered = [
        post_prob.index(post_prob_sorted[0]),
        round(post_prob_sorted[0], 3),
        post_prob.index(post_prob_sorted[1]),
        round(post_prob_sorted[1], 3),
    ]
    return cn_prob_filtered


def process_raw_call_gc(cn_prob, post_cutoff, keep_none=True):
    """
    Filter raw CN calls based on posterior probability cutoff.
    For gene conversion cases (SNVs between OTOA and pseudogene).
    
    Args:
        cn_prob: List of raw CN call results from call_reg1_cn
        post_cutoff: Posterior probability cutoff for confident calls
        keep_none: Whether to include None values in output
        
    Returns:
        List of filtered CN values
    """
    cn_prob_filtered = []
    for cn_call in cn_prob:
        if len(cn_call) == 1:
            call_value = cn_call[0]
            if call_value is not None or keep_none:
                cn_prob_filtered.append(call_value)
        
        # Resolve ambiguous CN=2 vs CN=3: prefer CN=2 when confidence is marginal
        elif (len(cn_call) > 2 and cn_call[0] == 3 and cn_call[2] == 2 
              and cn_call[1] < AMBIGUOUS_CALL_THRESHOLD):
            cn_prob_filtered.append(2)
        
        elif cn_call[1] > post_cutoff:
            cn_prob_filtered.append(cn_call[0])
        elif keep_none:
            cn_prob_filtered.append(None)
    
    return cn_prob_filtered


def process_raw_call_denovo(
    cn_prob, post_cutoff1, post_cutoff2, list_total_cn=None, keep_none=True
):
    """
    Filter raw CN calls based on posterior probability cutoff.
    For de novo variant calling (non-gene-conversion cases).
    
    For less confident calls that are not copy number zero,
    return the smaller CN call. Also keep the variant if called CN 
    equals total CN at the site (ensures variants are kept if CN>=1
    but we can't distinguish 1 vs 2).
    
    Args:
        cn_prob: List of raw CN call results
        post_cutoff1: Primary posterior probability cutoff
        post_cutoff2: Secondary cutoff for less confident calls
        list_total_cn: Optional list of total CN at each site
        keep_none: Whether to include None values
        
    Returns:
        List of filtered CN values
    """
    cn_prob_filtered = []
    for i, cn_call in enumerate(cn_prob):
        if len(cn_call) == 1:
            call_value = cn_call[0]
            if call_value is not None or keep_none:
                cn_prob_filtered.append(call_value)
        else:
            if list_total_cn is not None:
                total_cn = list_total_cn[i]
                keep_var = (cn_call[0] > 0 and cn_call[2] > 0) or (
                    cn_call[0] == total_cn or cn_call[2] == total_cn
                )
            else:
                keep_var = cn_call[0] > 0 and cn_call[2] > 0
            
            if cn_call[1] > post_cutoff1:
                cn_prob_filtered.append(cn_call[0])
            elif keep_var:
                if cn_call[1] > post_cutoff2:
                    cn_prob_filtered.append(cn_call[0])
                else:
                    cn_prob_filtered.append(min(cn_call[0], cn_call[2]))
            elif keep_none:
                cn_prob_filtered.append(None)
    
    return cn_prob_filtered
