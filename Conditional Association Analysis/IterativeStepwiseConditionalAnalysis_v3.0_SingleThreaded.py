"""
Iterative Stepwise Conditional Analysis using GCTA-COJO (SINGLE-THREADED VERSION)

This program performs iterative conditional analysis where:
1. Starts with a lead SNP
2. Runs conditional analysis
3. Identifies significant SNPs in a defined region
4. Adds them to the conditioning list
5. Repeats until no new significant SNPs are found

FEATURES:
- Single-threaded sequential processing (simpler, more stable)
- Automatic format detection for summary statistics
- Chromosome-specific COJO files for each lead SNP (uses CHR column from source)
- Works with both rsIDs and CHR:BP format SNP IDs
- Convergence visualization plots

Author: Based on ConditionalAnalysis.py
Date: October 2025
Version: 3.0 (Single-threaded, simplified)
"""

import os
import pandas as pd
import subprocess
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for cluster compatibility
from pathlib import Path
import time
import numpy as np


# ============================================================================
# SECTION 1: FORMAT DETECTION AND CONVERSION FUNCTIONS
# ============================================================================

def detect_file_format(summary_file):
    """
    Automatically detect whether the summary statistics file is from 
    PLINK logistic regression or meta-analysis.

    Detection logic:
    - Logistic files have: NMISS, TEST, STAT columns
    - Meta files have: P(R), OR(R), Q, I columns

    Parameters:
    -----------
    summary_file : str
        Path to summary statistics file

    Returns:
    --------
    str
        'logistic' or 'meta' or 'unknown'
    """
    try:
        # Read first few rows to check columns
        df = pd.read_csv(summary_file, nrows=10)
        columns = set(df.columns)

        # Check for logistic-specific columns
        logistic_markers = {'NMISS', 'TEST', 'STAT'}
        if logistic_markers.issubset(columns):
            print("Detected format: PLINK Logistic Regression")
            return 'logistic'

        # Check for meta-analysis-specific columns
        meta_markers = {'P(R)', 'Q', 'I'}
        if meta_markers.issubset(columns):
            print("Detected format: PLINK Meta-Analysis")
            return 'meta'

        # If COJO format already (SNP, A1, A2, freq, b, se, p, N)
        cojo_markers = {'SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N'}
        if cojo_markers.issubset(columns):
            print("Detected format: GCTA-COJO format (already converted)")
            return 'cojo'

        print("WARNING: Could not automatically detect file format")
        return 'unknown'

    except Exception as e:
        print(f"ERROR detecting file format: {str(e)}")
        return 'unknown'


def convert_to_cojo_logistic(logistic_file, output_file):
    """
    Convert PLINK logistic regression output to GCTA-COJO format.

    Required COJO format columns: SNP, A1, A2, freq, b, se, p, N
    PRESERVES: CHR, BP columns from source file

    Parameters:
    -----------
    logistic_file : str
        Path to PLINK logistic regression output file
    output_file : str
        Path for output COJO format file

    Returns:
    --------
    pd.DataFrame
        DataFrame in COJO format with CHR and BP columns
    """
    # Read the logistic file
    df = pd.read_csv(logistic_file)
    df.dropna(inplace=True)

    # Validate required columns exist
    required_cols = ['CHR', 'SNP', 'BP', 'A1_MINORALLELE', 'BETA', 'SE', 'P', 'NMISS']
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"Missing required column in logistic file: {col}")

    # Create new dataframe with required columns
    cojo_df = pd.DataFrame()
    cojo_df['CHR'] = df['CHR']  # PRESERVE chromosome column
    cojo_df['BP'] = df['BP']    # PRESERVE base pair column
    cojo_df['SNP'] = df['SNP']
    cojo_df['A1'] = df['A1_MINORALLELE']
    cojo_df['A2'] = df['A2_MAJORALLELE']
    cojo_df['freq'] = df['A1_MINORALLELEFREQ']
    cojo_df['b'] = df['BETA']
    cojo_df['se'] = df['SE']
    cojo_df['p'] = df['P']
    cojo_df['N'] = df['NMISS']

    # Save to tab-separated file (with CHR and BP for internal use)
    cojo_df.to_csv(output_file, sep='\t', index=False)
    print(f"Converted logistic file to COJO format: {output_file}")
    print(f"  Preserved CHR and BP columns from source file")

    return cojo_df


def convert_to_cojo_meta(meta_file, output_file, total_n):
    """
    Convert PLINK meta-analysis output to GCTA-COJO format.

    Required COJO format columns: SNP, A1, A2, freq, b, se, p, N
    PRESERVES: CHR, BP columns from source file

    Parameters:
    -----------
    meta_file : str
        Path to PLINK meta-analysis output file
    output_file : str
        Path for output COJO format file
    total_n : int
        Total sample size across all cohorts

    Returns:
    --------
    pd.DataFrame
        DataFrame in COJO format with CHR and BP columns
    """
    # Read the meta-analysis file
    df = pd.read_csv(meta_file)
    df.dropna(inplace=True)

    # Validate required columns exist
    required_cols = ['CHR', 'BP', 'SNP', 'A1_MinorAllele', 'A2_MajorAllele', 'BETA', 'SE', 'P']
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"Missing required column in meta file: {col}")

    # Create new dataframe with required columns
    cojo_df = pd.DataFrame()
    cojo_df['CHR'] = df['CHR']  # PRESERVE chromosome column
    cojo_df['BP'] = df['BP']    # PRESERVE base pair column
    cojo_df['SNP'] = df['SNP']
    cojo_df['A1'] = df['A1_MinorAllele']
    cojo_df['A2'] = df['A2_MajorAllele']
    cojo_df['freq'] = df['A1_MinorAlleleFreq']
    cojo_df['b'] = df['BETA']
    cojo_df['se'] = df['SE']
    cojo_df['p'] = df['P']
    cojo_df['N'] = total_n  # Use provided total sample size

    # Save to tab-separated file (with CHR and BP for internal use)
    cojo_df.to_csv(output_file, sep='\t', index=False)
    print(f"Converted meta-analysis file to COJO format: {output_file}")
    print(f"  Preserved CHR and BP columns from source file")

    return cojo_df


def auto_convert_to_cojo(summary_file, output_file, total_n=None):
    """
    Automatically detect format and convert to COJO format if needed.

    Parameters:
    -----------
    summary_file : str
        Path to summary statistics file
    output_file : str
        Path for output COJO format file
    total_n : int or None
        Total sample size (required for meta-analysis, optional for logistic)

    Returns:
    --------
    tuple
        (success: bool, cojo_df: pd.DataFrame or None)
    """
    file_format = detect_file_format(summary_file)

    if file_format == 'cojo':
        # Already in COJO format, check if it has CHR and BP columns
        df = pd.read_csv(summary_file, sep='\t')

        # Verify CHR and BP columns exist
        if 'CHR' not in df.columns or 'BP' not in df.columns:
            print("ERROR: COJO file missing CHR and/or BP columns")
            return False, None

        df.to_csv(output_file, sep='\t', index=False)
        print(f"File already in COJO format with CHR and BP columns. Copied to: {output_file}")
        return True, df

    elif file_format == 'logistic':
        cojo_df = convert_to_cojo_logistic(summary_file, output_file)
        return True, cojo_df

    elif file_format == 'meta':
        if total_n is None:
            print("ERROR: Meta-analysis format detected but total_n not provided")
            return False, None
        cojo_df = convert_to_cojo_meta(summary_file, output_file, total_n)
        return True, cojo_df

    else:
        print("ERROR: Unknown file format. Manual conversion required.")
        return False, None


def create_chromosome_specific_cojo(cojo_df, chromosome, output_file):
    """
    Create a chromosome-specific COJO file containing only SNPs from specified chromosome.

    Uses the CHR column from the source file (not parsing SNP IDs).
    This works with both rsIDs (rs1234) and CHR:BP format (6:31768055).

    Parameters:
    -----------
    cojo_df : pd.DataFrame
        Full COJO format dataframe with CHR column
    chromosome : int
        Chromosome number to filter
    output_file : str
        Path for output chromosome-specific COJO file

    Returns:
    --------
    pd.DataFrame
        Filtered DataFrame containing only SNPs from specified chromosome
    """
    # Verify CHR column exists
    if 'CHR' not in cojo_df.columns:
        raise ValueError("COJO dataframe missing CHR column. Cannot filter by chromosome.")

    # Filter to specified chromosome using CHR column from source file
    chr_df = cojo_df[cojo_df['CHR'] == chromosome].copy()

    # Create output dataframe WITHOUT CHR and BP columns (GCTA doesn't need them)
    # GCTA only needs: SNP, A1, A2, freq, b, se, p, N
    chr_df_output = chr_df[['SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N']].copy()

    # Save to tab-separated file
    chr_df_output.to_csv(output_file, sep='\t', index=False)

    print(f"Created chromosome {chromosome}-specific COJO file with {len(chr_df)} SNPs: {output_file}")
    print(f"  Filtered using CHR column from source file (works with rsIDs)")

    return chr_df


# ============================================================================
# SECTION 2: SNP LIST MANAGEMENT FUNCTIONS
# ============================================================================

def create_snp_list_file(snp_list, output_file):
    """
    Create a text file containing SNP IDs for GCTA --cojo-cond option.

    Parameters:
    -----------
    snp_list : list
        List of SNP IDs to condition on
    output_file : str
        Path for output SNP list file

    Returns:
    --------
    None
    """
    with open(output_file, 'w') as f:
        for snp in snp_list:
            f.write(f"{snp}\n")
    print(f"Created SNP list file with {len(snp_list)} SNPs: {output_file}")


def read_lead_snps_file(lead_snps_file):
    """
    Read the CSV file containing starting lead SNPs.

    Expected columns: CHR, SNP, BP

    Parameters:
    -----------
    lead_snps_file : str
        Path to CSV file with lead SNPs

    Returns:
    --------
    pd.DataFrame
        DataFrame with CHR, SNP, BP columns
    """
    df = pd.read_csv(lead_snps_file)

    # Validate required columns
    required_cols = ['CHR', 'SNP', 'BP']
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"Missing required column in lead SNPs file: {col}")

    print(f"Read {len(df)} lead SNPs from {lead_snps_file}")
    return df


# ============================================================================
# SECTION 3: GCTA COMMAND EXECUTION FUNCTIONS
# ============================================================================

def build_gcta_command(gcta_path, bfile, cojo_file, snp_list_file, 
                       chromosome, output_prefix, maf, window, 
                       collinearity, diff_freq, thread_num, use_gc):
    """
    Build GCTA-COJO conditional analysis command.

    Parameters:
    -----------
    gcta_path : str
        Path to GCTA executable
    bfile : str
        Prefix for PLINK binary files (.bed/.bim/.fam)
    cojo_file : str
        Path to COJO format summary statistics file
    snp_list_file : str
        Path to file containing SNPs to condition on
    chromosome : int
        Chromosome number
    output_prefix : str
        Prefix for output files
    maf : float
        Minor allele frequency threshold
    window : int
        LD window size in kb
    collinearity : float
        Collinearity threshold
    diff_freq : float
        Allele frequency difference threshold
    thread_num : int
        Number of threads to use
    use_gc : bool
        Whether to use genomic control correction

    Returns:
    --------
    str
        Complete GCTA command
    """
    cmd = f"{gcta_path} --bfile {bfile} --thread-num {thread_num} "
    cmd += f"--cojo-file {cojo_file} --maf {maf} --cojo-wind {window} "
    cmd += f"--cojo-collinear {collinearity} --diff-freq {diff_freq} "
    cmd += f"--chr {chromosome} --out {output_prefix} "
    cmd += f"--cojo-cond {snp_list_file}"

    if use_gc:
        cmd += " --cojo-gc"

    return cmd


def execute_gcta_command(command):
    """
    Execute a single GCTA command using subprocess.

    Parameters:
    -----------
    command : str
        Complete GCTA command to execute

    Returns:
    --------
    int
        Return code (0 for success, non-zero for failure)
    """
    try:
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"ERROR: Command failed with return code {result.returncode}")
            print(f"STDERR: {result.stderr}")
        return result.returncode
    except Exception as e:
        print(f"ERROR executing command: {str(e)}")
        return -1


# ============================================================================
# SECTION 4: RESULTS PARSING AND FILTERING FUNCTIONS
# ============================================================================

def parse_cojo_results(cojo_output_file):
    """
    Parse GCTA-COJO conditional analysis results file (.cma.cojo).

    Expected columns: Chr, SNP, bp, refA, freq, b, se, p, n, freq_geno, bC, bC_se, pC

    Parameters:
    -----------
    cojo_output_file : str
        Path to .cma.cojo output file

    Returns:
    --------
    pd.DataFrame or None
        DataFrame with results, or None if file doesn't exist
    """
    if not os.path.exists(cojo_output_file):
        print(f"WARNING: Results file not found: {cojo_output_file}")
        return None

    try:
        # GCTA outputs space-separated files
        df = pd.read_csv(cojo_output_file, sep='\s+')
        print(f"Parsed {len(df)} SNPs from {cojo_output_file}")
        return df
    except Exception as e:
        print(f"ERROR parsing results file: {str(e)}")
        return None


def filter_snps_by_region(results_df, lead_chr, lead_bp, window_kb):
    """
    Filter SNPs to those on the same chromosome and within ±window_kb of lead SNP.

    Parameters:
    -----------
    results_df : pd.DataFrame
        Conditional analysis results
    lead_chr : int
        Chromosome of lead SNP
    lead_bp : int
        Base pair position of lead SNP
    window_kb : int
        Window size in kilobases (e.g., 500 for ±500kb)

    Returns:
    --------
    pd.DataFrame
        Filtered DataFrame
    """
    # Convert window from kb to bp
    window_bp = window_kb * 1000

    # Calculate region boundaries
    region_start = lead_bp - window_bp
    region_end = lead_bp + window_bp

    # Filter by chromosome and position
    # GCTA output has 'Chr' and 'bp' columns
    filtered = results_df[
        (results_df['Chr'] == lead_chr) & 
        (results_df['bp'] >= region_start) & 
        (results_df['bp'] <= region_end)
    ].copy()

    print(f"Filtered to {len(filtered)} SNPs in chr{lead_chr}:{region_start}-{region_end}")
    return filtered


def extract_significant_snps(results_df, p_threshold):
    """
    Extract SNPs with conditional p-value (pC) below threshold.

    Parameters:
    -----------
    results_df : pd.DataFrame
        Conditional analysis results
    p_threshold : float
        P-value threshold (e.g., 5e-6)

    Returns:
    --------
    list
        List of significant SNP IDs
    """
    significant = results_df[results_df['pC'] <= p_threshold]
    snp_list = significant['SNP'].tolist()

    print(f"Found {len(snp_list)} SNPs with pC <= {p_threshold}")
    return snp_list


def sort_and_save_results(results_df, output_file, sort_by='pC'):
    """
    Sort results by specified column and save to CSV.

    Parameters:
    -----------
    results_df : pd.DataFrame
        Conditional analysis results
    output_file : str
        Path for output CSV file
    sort_by : str
        Column name to sort by (default: 'pC')

    Returns:
    --------
    None
    """
    sorted_df = results_df.sort_values(by=sort_by)
    sorted_df.to_csv(output_file, index=False)
    print(f"Saved sorted results to {output_file}")


# ============================================================================
# SECTION 5: VISUALIZATION FUNCTIONS
# ============================================================================

def create_convergence_plot(iteration_history, output_dir, lead_snp):
    """
    Create visualization showing convergence of iterative analysis.

    Creates two plots:
    1. Number of conditioning SNPs per iteration
    2. Number of new significant SNPs found per iteration

    Parameters:
    -----------
    iteration_history : list of dict
        History of each iteration
    output_dir : str
        Directory for output files
    lead_snp : str
        Lead SNP identifier for plot title

    Returns:
    --------
    None
    """
    # Extract data from history
    iterations = [h['iteration'] for h in iteration_history]
    conditioning_snps = [h['num_conditioning_snps'] for h in iteration_history]
    new_snps = [h['new_snps_found'] for h in iteration_history]

    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

    # Plot 1: Cumulative conditioning SNPs
    ax1.plot(iterations, conditioning_snps, marker='o', linewidth=2, markersize=8, color='steelblue')
    ax1.set_xlabel('Iteration', fontsize=12)
    ax1.set_ylabel('Number of Conditioning SNPs', fontsize=12)
    ax1.set_title(f'Convergence Analysis: {lead_snp}\nCumulative Conditioning SNPs', 
                  fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.set_xticks(iterations)

    # Plot 2: New SNPs found per iteration
    ax2.bar(iterations, new_snps, color='coral', alpha=0.7)
    ax2.set_xlabel('Iteration', fontsize=12)
    ax2.set_ylabel('New Significant SNPs Found', fontsize=12)
    ax2.set_title('New Significant SNPs per Iteration', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3, axis='y')
    ax2.set_xticks(iterations)

    # Add horizontal line at y=0 for second plot
    ax2.axhline(y=0, color='red', linestyle='--', linewidth=1, alpha=0.5)

    plt.tight_layout()

    # Save plot
    plot_file = os.path.join(output_dir, f'Convergence_Plot_{lead_snp}.png')
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Convergence plot saved to {plot_file}")


def create_pvalue_progression_plot(iteration_results_files, output_dir, lead_snp):
    """
    Create plot showing how p-values change across iterations for top SNPs.

    Parameters:
    -----------
    iteration_results_files : list
        List of tuples (iteration_num, results_file_path)
    output_dir : str
        Directory for output files
    lead_snp : str
        Lead SNP identifier

    Returns:
    --------
    None
    """
    # Track p-values for top SNPs across iterations
    snp_pvalues = {}
    iterations = []

    for iter_num, results_file in iteration_results_files:
        df = parse_cojo_results(results_file)
        if df is None:
            continue

        iterations.append(iter_num)

        # Get top 10 SNPs by pC
        top_snps = df.nsmallest(10, 'pC')

        for _, row in top_snps.iterrows():
            snp_id = row['SNP']
            pC = row['pC']

            if snp_id not in snp_pvalues:
                snp_pvalues[snp_id] = {}
            snp_pvalues[snp_id][iter_num] = pC

    if len(snp_pvalues) == 0:
        print("No SNPs to plot for p-value progression")
        return

    # Create plot
    fig, ax = plt.subplots(figsize=(12, 8))

    for snp_id, pvals in snp_pvalues.items():
        iters = sorted(pvals.keys())
        pval_list = [pvals[i] for i in iters]

        # Plot on log scale
        ax.plot(iters, [-np.log10(p) for p in pval_list], 
                marker='o', label=snp_id, alpha=0.7)

    ax.set_xlabel('Iteration', fontsize=12)
    ax.set_ylabel('-log10(pC)', fontsize=12)
    ax.set_title(f'P-value Progression Across Iterations: {lead_snp}', 
                 fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)

    plt.tight_layout()

    # Save plot
    plot_file = os.path.join(output_dir, f'PValue_Progression_{lead_snp}.png')
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"P-value progression plot saved to {plot_file}")


# ============================================================================
# SECTION 6: ITERATION MANAGEMENT FUNCTIONS
# ============================================================================

def run_conditional_iteration(gcta_path, bfile, cojo_file, snp_list, 
                              lead_chr, output_dir, iteration_num, 
                              maf, window, collinearity, diff_freq, 
                              thread_num, use_gc):
    """
    Run a single iteration of conditional analysis.

    Parameters:
    -----------
    gcta_path : str
        Path to GCTA executable
    bfile : str
        Prefix for PLINK binary files
    cojo_file : str
        Path to COJO format summary statistics
    snp_list : list
        List of SNPs to condition on
    lead_chr : int
        Chromosome number
    output_dir : str
        Directory for output files
    iteration_num : int
        Current iteration number
    maf : float
        Minor allele frequency threshold
    window : int
        LD window size in kb
    collinearity : float
        Collinearity threshold
    diff_freq : float
        Allele frequency difference threshold
    thread_num : int
        Number of threads
    use_gc : bool
        Whether to use genomic control

    Returns:
    --------
    tuple
        (success: bool, results_file: str)
    """
    # Create SNP list file for this iteration
    snp_list_file = os.path.join(output_dir, f"Iteration_{iteration_num}_ConditionedSNPs.txt")
    create_snp_list_file(snp_list, snp_list_file)

    # Build output prefix
    output_prefix = os.path.join(output_dir, f"Iteration_{iteration_num}_Conditional")

    # Build and execute GCTA command
    command = build_gcta_command(
        gcta_path, bfile, cojo_file, snp_list_file, lead_chr, 
        output_prefix, maf, window, collinearity, diff_freq, 
        thread_num, use_gc
    )

    print(f"\n{'='*70}")
    print(f"Running Iteration {iteration_num} with {len(snp_list)} conditioning SNPs")
    print(f"{'='*70}")
    print(f"Command: {command}\n")

    # Execute command
    return_code = execute_gcta_command(command)

    # Check if successful
    results_file = output_prefix + ".cma.cojo"
    success = (return_code == 0) and os.path.exists(results_file)

    if success:
        print(f"Iteration {iteration_num} completed successfully")
    else:
        print(f"ERROR: Iteration {iteration_num} failed")

    return success, results_file


def process_iteration_results(results_file, output_dir, iteration_num, 
                              lead_chr, lead_bp, window_kb, p_threshold):
    """
    Process results from a single iteration: sort, filter, and extract significant SNPs.

    Parameters:
    -----------
    results_file : str
        Path to .cma.cojo results file
    output_dir : str
        Directory for output files
    iteration_num : int
        Current iteration number
    lead_chr : int
        Chromosome of lead SNP
    lead_bp : int
        Base pair position of lead SNP
    window_kb : int
        Window size in kilobases
    p_threshold : float
        P-value threshold for significance

    Returns:
    --------
    tuple
        (new_snps: list, sorted_results_file: str)
    """
    # Parse results
    results_df = parse_cojo_results(results_file)
    if results_df is None:
        return [], None

    # Sort and save full results
    sorted_file = os.path.join(output_dir, f"Iteration_{iteration_num}_Conditional_Sorted.csv")
    sort_and_save_results(results_df, sorted_file)

    # Filter to region of interest
    filtered_df = filter_snps_by_region(results_df, lead_chr, lead_bp, window_kb)

    # Save filtered results
    filtered_file = os.path.join(output_dir, f"Iteration_{iteration_num}_RegionOfInterest.csv")
    filtered_df.to_csv(filtered_file, index=False)
    print(f"Saved region-of-interest results to {filtered_file}")

    # Extract significant SNPs
    new_snps = extract_significant_snps(filtered_df, p_threshold)

    return new_snps, sorted_file


def create_iteration_summary(output_dir, lead_snp, lead_chr, lead_bp, 
                            total_iterations, iteration_history, 
                            final_snp_list, p_threshold, window_kb):
    """
    Create a summary file of the iterative analysis.

    Parameters:
    -----------
    output_dir : str
        Directory for output files
    lead_snp : str
        Lead SNP identifier
    lead_chr : int
        Chromosome of lead SNP
    lead_bp : int
        Base pair position of lead SNP
    total_iterations : int
        Total number of iterations performed
    iteration_history : list of dict
        History of each iteration with details
    final_snp_list : list
        Final list of conditioning SNPs
    p_threshold : float
        P-value threshold used
    window_kb : int
        Window size used

    Returns:
    --------
    None
    """
    summary_file = os.path.join(output_dir, "Analysis_Summary.txt")

    with open(summary_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("ITERATIVE STEPWISE CONDITIONAL ANALYSIS SUMMARY\n")
        f.write("="*80 + "\n\n")

        f.write(f"Lead SNP: {lead_snp}\n")
        f.write(f"Chromosome: {lead_chr}\n")
        f.write(f"Base Pair Position: {lead_bp:,}\n")
        f.write(f"Region of Interest: ±{window_kb} kb\n")
        f.write(f"P-value Threshold: {p_threshold}\n")
        f.write(f"Total Iterations: {total_iterations}\n\n")

        f.write("="*80 + "\n")
        f.write("ITERATION HISTORY\n")
        f.write("="*80 + "\n\n")

        for i, hist in enumerate(iteration_history, 1):
            f.write(f"Iteration {i}:\n")
            f.write(f"  Conditioning SNPs: {hist['num_conditioning_snps']}\n")
            f.write(f"  New Significant SNPs Found: {hist['new_snps_found']}\n")
            if hist['new_snps_found'] > 0:
                f.write(f"  New SNPs: {', '.join(hist['new_snps'])}\n")
            f.write("\n")

        f.write("="*80 + "\n")
        f.write("FINAL RESULTS\n")
        f.write("="*80 + "\n\n")

        f.write(f"Total Conditioning SNPs: {len(final_snp_list)}\n")
        f.write(f"\nFinal SNP List:\n")
        for snp in final_snp_list:
            f.write(f"  - {snp}\n")

    print(f"\nSummary saved to {summary_file}")


# ============================================================================
# SECTION 7: MAIN ITERATIVE ANALYSIS FUNCTION
# ============================================================================

def run_iterative_conditional_analysis(gcta_path, bfile, master_cojo_df, 
                                      lead_snp, lead_chr, lead_bp, 
                                      output_base_dir, p_threshold, 
                                      window_kb, maf, ld_window, 
                                      collinearity, diff_freq, 
                                      thread_num, use_gc, create_plots=True):
    """
    Run complete iterative stepwise conditional analysis for a single lead SNP.

    This function performs the following steps:
    1. Create chromosome-specific COJO file for this lead SNP (using CHR column)
    2. Start with the lead SNP
    3. Run conditional analysis
    4. Identify significant SNPs in the region of interest
    5. Add them to the conditioning list
    6. Repeat until no new significant SNPs are found
    7. Create visualization plots

    Parameters:
    -----------
    gcta_path : str
        Path to GCTA executable
    bfile : str
        Prefix for PLINK binary files
    master_cojo_df : pd.DataFrame
        Full COJO format dataframe (all chromosomes) with CHR column
    lead_snp : str
        Lead SNP identifier
    lead_chr : int
        Chromosome of lead SNP
    lead_bp : int
        Base pair position of lead SNP
    output_base_dir : str
        Base directory for all outputs
    p_threshold : float
        P-value threshold for significance (e.g., 5e-6)
    window_kb : int
        Window size in kilobases for region of interest (e.g., 500)
    maf : float
        Minor allele frequency threshold
    ld_window : int
        LD window size in kb for GCTA
    collinearity : float
        Collinearity threshold for GCTA
    diff_freq : float
        Allele frequency difference threshold for GCTA
    thread_num : int
        Number of threads to use
    use_gc : bool
        Whether to use genomic control correction
    create_plots : bool
        Whether to create visualization plots

    Returns:
    --------
    dict
        Dictionary with analysis results and statistics
    """
    # Create output directory for this lead SNP
    output_dir = os.path.join(output_base_dir, lead_snp)
    os.makedirs(output_dir, exist_ok=True)
    print(f"\n{'#'*80}")
    print(f"# Starting iterative analysis for lead SNP: {lead_snp}")
    print(f"# Chromosome: {lead_chr}")
    print(f"# Output directory: {output_dir}")
    print(f"{'#'*80}\n")

    # Create chromosome-specific COJO file for this lead SNP
    # Uses CHR column from source file - works with both rsIDs and CHR:BP format
    cojo_file = os.path.join(output_dir, f"{lead_snp}_Chr{lead_chr}_COJO.txt")

    try:
        chr_cojo_df = create_chromosome_specific_cojo(master_cojo_df, lead_chr, cojo_file)
    except Exception as e:
        print(f"ERROR creating chromosome-specific COJO file: {str(e)}")
        return {
            'lead_snp': lead_snp,
            'total_iterations': 0,
            'final_snp_count': 0,
            'final_snp_list': [],
            'success': False,
            'error': f'Failed to create chromosome-specific COJO file: {str(e)}'
        }

    if len(chr_cojo_df) == 0:
        print(f"ERROR: No SNPs found on chromosome {lead_chr} in COJO file")
        return {
            'lead_snp': lead_snp,
            'total_iterations': 0,
            'final_snp_count': 0,
            'final_snp_list': [],
            'success': False,
            'error': f'No SNPs on chromosome {lead_chr}'
        }

    # Initialize
    current_snp_list = [lead_snp]
    iteration_num = 1
    iteration_history = []
    iteration_results_files = []

    # Main iteration loop
    while True:
        # Run conditional analysis for this iteration
        success, results_file = run_conditional_iteration(
            gcta_path, bfile, cojo_file, current_snp_list, lead_chr, 
            output_dir, iteration_num, maf, ld_window, collinearity, 
            diff_freq, thread_num, use_gc
        )

        if not success:
            print(f"\nERROR: Iteration {iteration_num} failed. Stopping analysis.")
            break

        # Store results file for plotting
        iteration_results_files.append((iteration_num, results_file))

        # Process results
        new_snps, sorted_file = process_iteration_results(
            results_file, output_dir, iteration_num, lead_chr, 
            lead_bp, window_kb, p_threshold
        )

        # Remove SNPs that are already in the conditioning list
        new_snps = [snp for snp in new_snps if snp not in current_snp_list]

        # Record iteration history
        iteration_history.append({
            'iteration': iteration_num,
            'num_conditioning_snps': len(current_snp_list),
            'new_snps_found': len(new_snps),
            'new_snps': new_snps
        })

        print(f"\nIteration {iteration_num} Summary:")
        print(f"  - Conditioned on {len(current_snp_list)} SNPs")
        print(f"  - Found {len(new_snps)} new significant SNPs")

        # Check convergence
        if len(new_snps) == 0:
            print(f"\n{'='*70}")
            print(f"CONVERGENCE REACHED: No new significant SNPs found")
            print(f"{'='*70}\n")
            break

        # Add new SNPs to conditioning list
        current_snp_list.extend(new_snps)
        print(f"  - Total conditioning SNPs now: {len(current_snp_list)}\n")

        # Increment iteration counter
        iteration_num += 1

    # Create summary
    create_iteration_summary(
        output_dir, lead_snp, lead_chr, lead_bp, iteration_num, 
        iteration_history, current_snp_list, p_threshold, window_kb
    )

    # Create visualization plots
    if create_plots and len(iteration_history) > 0:
        try:
            create_convergence_plot(iteration_history, output_dir, lead_snp)
            if len(iteration_results_files) > 1:
                create_pvalue_progression_plot(iteration_results_files, output_dir, lead_snp)
        except Exception as e:
            print(f"WARNING: Could not create plots: {str(e)}")

    return {
        'lead_snp': lead_snp,
        'total_iterations': iteration_num,
        'final_snp_count': len(current_snp_list),
        'final_snp_list': current_snp_list,
        'success': True
    }


# ============================================================================
# SECTION 8: SEQUENTIAL BATCH PROCESSING FUNCTION
# ============================================================================

def sequential_process_lead_snps(gcta_path, bfile, master_cojo_df, lead_snps_df, 
                                 output_base_dir, p_threshold, window_kb, 
                                 maf, ld_window, collinearity, diff_freq, 
                                 thread_num, use_gc, create_plots=True):
    """
    Process multiple lead SNPs sequentially (one after another).

    Each lead SNP is processed completely before moving to the next one.
    Much simpler than parallel processing - no threading complications.

    Parameters:
    -----------
    gcta_path : str
        Path to GCTA executable
    bfile : str
        Prefix for PLINK binary files
    master_cojo_df : pd.DataFrame
        Full COJO format dataframe (all chromosomes) with CHR column
    lead_snps_df : pd.DataFrame
        DataFrame with CHR, SNP, BP columns for lead SNPs
    output_base_dir : str
        Base directory for all outputs
    p_threshold : float
        P-value threshold for significance
    window_kb : int
        Window size in kilobases for region of interest
    maf : float
        Minor allele frequency threshold
    ld_window : int
        LD window size in kb for GCTA
    collinearity : float
        Collinearity threshold for GCTA
    diff_freq : float
        Allele frequency difference threshold for GCTA
    thread_num : int
        Number of threads for GCTA to use
    use_gc : bool
        Whether to use genomic control correction
    create_plots : bool
        Whether to create visualization plots

    Returns:
    --------
    list of dict
        List of results for each lead SNP
    """
    total_snps = len(lead_snps_df)
    results_list = []

    print(f"\n{'#'*80}")
    print(f"# SEQUENTIAL PROCESSING: {total_snps} LEAD SNPs")
    print(f"# Each lead SNP will be processed one after another")
    print(f"# GCTA threads per SNP: {thread_num}")
    print(f"{'#'*80}\n")

    start_time = time.time()

    # Process each lead SNP one by one
    for idx, row in lead_snps_df.iterrows():
        lead_snp = row['SNP']
        lead_chr = int(row['CHR'])
        lead_bp = int(row['BP'])

        print(f"\n{'='*80}")
        print(f"[{idx+1}/{total_snps}] Processing lead SNP: {lead_snp}")
        print(f"{'='*80}")

        try:
            result = run_iterative_conditional_analysis(
                gcta_path, bfile, master_cojo_df, lead_snp, lead_chr, lead_bp, 
                output_base_dir, p_threshold, window_kb, maf, ld_window, 
                collinearity, diff_freq, thread_num, use_gc, create_plots
            )
            results_list.append(result)

            print(f"\n[{idx+1}/{total_snps}] COMPLETED: {lead_snp}")
            print(f"  Iterations: {result['total_iterations']}")
            print(f"  Final SNPs: {result['final_snp_count']}")

        except Exception as e:
            print(f"\nERROR processing {lead_snp}: {str(e)}")
            import traceback
            traceback.print_exc()
            results_list.append({
                'lead_snp': lead_snp,
                'total_iterations': 0,
                'final_snp_count': 0,
                'final_snp_list': [],
                'success': False,
                'error': str(e)
            })

    elapsed_time = time.time() - start_time

    print(f"\n{'#'*80}")
    print(f"# ALL LEAD SNPs PROCESSED")
    print(f"# Total time: {elapsed_time/60:.2f} minutes")
    print(f"# Average time per SNP: {elapsed_time/total_snps:.1f} seconds")
    print(f"{'#'*80}\n")

    # Create batch summary
    create_batch_summary(results_list, output_base_dir)

    return results_list


def create_batch_summary(results_list, output_base_dir):
    """
    Create a summary file for batch processing of multiple lead SNPs.

    Parameters:
    -----------
    results_list : list of dict
        List of results for each lead SNP
    output_base_dir : str
        Base directory for outputs

    Returns:
    --------
    None
    """
    summary_file = os.path.join(output_base_dir, "Batch_Processing_Summary.txt")

    with open(summary_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("BATCH PROCESSING SUMMARY (SEQUENTIAL)\n")
        f.write("="*80 + "\n\n")

        f.write(f"Total Lead SNPs Processed: {len(results_list)}\n")
        successful = sum(1 for r in results_list if r['success'])
        failed = len(results_list) - successful
        f.write(f"Successful: {successful}\n")
        f.write(f"Failed: {failed}\n\n")

        f.write("="*80 + "\n")
        f.write("INDIVIDUAL SNP RESULTS\n")
        f.write("="*80 + "\n\n")

        for result in results_list:
            f.write(f"Lead SNP: {result['lead_snp']}\n")
            f.write(f"  Total Iterations: {result['total_iterations']}\n")
            f.write(f"  Final Conditioning SNPs: {result['final_snp_count']}\n")
            f.write(f"  Success: {result['success']}\n")
            if not result['success'] and 'error' in result:
                f.write(f"  Error: {result['error']}\n")
            f.write("\n")

    print(f"\nBatch summary saved to {summary_file}")


# ============================================================================
# SECTION 9: MAIN EXECUTION FUNCTION
# ============================================================================

def main():
    """
    Main execution function with user input collection.

    This function collects all required parameters from the user and 
    runs the iterative stepwise conditional analysis sequentially.
    """
    print("="*80)
    print("ITERATIVE STEPWISE CONDITIONAL ANALYSIS v3.0 (SINGLE-THREADED)")
    print("="*80)
    print()

    # Collect user inputs
    print("Please provide the following information:\n")

    # GCTA path
    gcta_path = input("Path to GCTA executable: ").strip()

    # Genotype files
    bfile = input("Path to PLINK binary file prefix (without .bed/.bim/.fam): ").strip()

    # Summary statistics with auto-detection
    summary_file = input("\nPath to summary statistics file: ").strip()

    # Create temporary master COJO file
    temp_master_cojo = summary_file.replace('.csv', '_Master_COJO_with_CHR_BP.TabSep')

    # Auto-detect format
    file_format = detect_file_format(summary_file)

    if file_format == 'meta':
        total_n = int(input("Total sample size across all cohorts: ").strip())
        success, master_cojo_df = auto_convert_to_cojo(summary_file, temp_master_cojo, total_n)
    elif file_format == 'unknown':
        print("\nCould not auto-detect format. Please specify:")
        print("  1. PLINK logistic regression output")
        print("  2. PLINK meta-analysis output")
        print("  3. Already in COJO format (must have CHR and BP columns)")
        file_type = input("Enter choice (1, 2, or 3): ").strip()

        if file_type == '1':
            master_cojo_df = convert_to_cojo_logistic(summary_file, temp_master_cojo)
            success = True
        elif file_type == '2':
            total_n = int(input("Total sample size across all cohorts: ").strip())
            master_cojo_df = convert_to_cojo_meta(summary_file, temp_master_cojo, total_n)
            success = True
        elif file_type == '3':
            master_cojo_df = pd.read_csv(summary_file, sep='\t')
            if 'CHR' not in master_cojo_df.columns or 'BP' not in master_cojo_df.columns:
                print("ERROR: File missing CHR and/or BP columns. Exiting.")
                return
            success = True
        else:
            print("Invalid choice. Exiting.")
            return
    else:
        success, master_cojo_df = auto_convert_to_cojo(summary_file, temp_master_cojo)

    if not success or master_cojo_df is None:
        print("ERROR: Could not convert summary statistics to COJO format. Exiting.")
        return

    print(f"\nMaster COJO file created with {len(master_cojo_df)} SNPs")
    print(f"  Contains CHR and BP columns from source file")

    # Lead SNPs file
    lead_snps_file = input("\nPath to lead SNPs CSV file (with CHR, SNP, BP columns): ").strip()
    lead_snps_df = read_lead_snps_file(lead_snps_file)

    # Output directory
    output_base_dir = input("\nOutput base directory: ").strip()
    os.makedirs(output_base_dir, exist_ok=True)

    # Analysis parameters
    print("\nAnalysis parameters:")
    p_threshold = float(input("P-value threshold (e.g., 5e-6): ").strip())
    window_kb = int(input("Region window size in kb (e.g., 500): ").strip())

    # GCTA parameters
    print("\nGCTA parameters (press Enter for defaults):")
    maf = input("MAF threshold [0.005]: ").strip()
    maf = float(maf) if maf else 0.005

    ld_window = input("LD window in kb [10000]: ").strip()
    ld_window = int(ld_window) if ld_window else 10000

    collinearity = input("Collinearity threshold [0.9]: ").strip()
    collinearity = float(collinearity) if collinearity else 0.9

    diff_freq = input("Allele frequency difference threshold [0.2]: ").strip()
    diff_freq = float(diff_freq) if diff_freq else 0.2

    thread_num = input("Number of threads for GCTA [20]: ").strip()
    thread_num = int(thread_num) if thread_num else 20

    # Genomic control option
    use_gc_input = input("\nUse genomic control correction? (yes/no) [no]: ").strip().lower()
    use_gc = use_gc_input in ['yes', 'y']

    # Visualization option
    create_plots_input = input("Create convergence plots? (yes/no) [yes]: ").strip().lower()
    create_plots = create_plots_input not in ['no', 'n']

    print("\n" + "="*80)
    print("Starting sequential analysis...")
    print("="*80 + "\n")

    # Run sequential processing
    results = sequential_process_lead_snps(
        gcta_path, bfile, master_cojo_df, lead_snps_df, output_base_dir, 
        p_threshold, window_kb, maf, ld_window, collinearity, 
        diff_freq, thread_num, use_gc, create_plots
    )

    print("\n" + "="*80)
    print("ANALYSIS COMPLETE")
    print("="*80)
    print(f"\nResults saved to: {output_base_dir}")
    print(f"Total lead SNPs processed: {len(results)}")
    successful = sum(1 for r in results if r['success'])
    print(f"Successful: {successful}")
    print(f"Failed: {len(results) - successful}")


if __name__ == "__main__":
    main()
