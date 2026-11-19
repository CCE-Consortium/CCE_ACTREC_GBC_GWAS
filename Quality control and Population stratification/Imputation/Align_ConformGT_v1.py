#!/usr/bin/env python3
"""
Conform-GT Alignment Pipeline - Simple Version
Converts PLINK BED to VCF once, then runs Conform-GT for all chromosomes in parallel batches.
Log processing is serial.
"""

import os
import subprocess
import threading
import pandas as pd
from pathlib import Path
import math

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def convert_plink_to_vcf(plink_prefix, output_vcf, plinkPath, plinkPath2):
    """Convert PLINK BED fileset to VCF (all chromosomes)."""
    os.system(f"{plinkPath} --bfile {plink_prefix} --list-duplicate-vars ids-only suppress-first --out {plink_prefix}_DUPS")
    os.system(f"{plinkPath} --bfile {plink_prefix} --exclude {plink_prefix}_DUPS.dupvar --make-bed --out {plink_prefix}_DupsCleaned")
    os.system(f"{plinkPath2} --bfile {plink_prefix}_DupsCleaned --export vcf bgz --out {output_vcf}")
    print(f"Converted PLINK to VCF: {output_vcf}")


def run_conform_gt(conform_gt_jar, reference_vcf, target_vcf, chromosome, output_prefix):
    """Run Conform-GT for a specific chromosome."""
    output_full = f"{output_prefix}_aligned"
    os.system(f"java -jar {conform_gt_jar} ref={reference_vcf} gt={target_vcf} chrom={chromosome} out={output_full} match=POS strict=false")
    # result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    # if result.returncode == 0:
    #     print(f"Chr {chromosome}: Conform-GT aligned")
    # else:
    #     print(f"Chr {chromosome}: Conform-GT failed")


def threads_list_run(ThreadsList):
    """Run all threads in list, wait for completion."""
    for t in ThreadsList:
        t.start()
    for t in ThreadsList:
        t.join()


def run_alignment_parallel(reference_vcf, target_vcf, conform_gt_jar, 
                           output_dir, output_prefix, max_workers):
    """Run Conform-GT alignment for all chromosomes in parallel batches."""
    
    print(f"Starting alignment with {max_workers} workers...")
    
    ThreadsList = []
    for chromosome in range(1, 27):
        output_path = os.path.join(output_dir, f"{output_prefix}_{chromosome}")
        
        ThreadObject = threading.Thread(
            target=run_conform_gt,
            args=(conform_gt_jar, reference_vcf, target_vcf, chromosome, output_path)
        )
        ThreadsList.append(ThreadObject)
        
        # # Run in batches of max_workers
        # if len(ThreadsList) <= max_workers:
        #     threads_list_run(ThreadsList)
        # else:
        #     for i in range(1,math.ceil(22/max_workers)):
        #         threads_list_run(ThreadsList[])
    # max_workers = 8     # this is where we are deciding the total number of alignment chromosomes that will be allowed to run simultaneuosly
    
    # algo for the below nested for loops..
    # find out how many batches will be run (number of chromosomes divided by the total number of worker threads that we can afford)
        #this gives us the  number of chromosomes that will be aligned in parallel
    # we need to get exactly that many threads out of the ThreadsList, and run it in one batch.
    # Therefore, for every batch....for X number of Chromosomes that we are going to align in one batch in parallel.....pop() those many threads from the Threads List and put it in a new batch list.
    # once we have multiple lists of batches.....we run those batches.
    num_of_chrs_in_batch = math.ceil(26/max_workers)
    batches = list()
    for i in range(0,max_workers):
        newBatch = list()
        for j in range(1, num_of_chrs_in_batch+1):
            try:
                newBatch.append(ThreadsList.pop())
            except IndexError:
                break
        batches.append(newBatch)
        
    for b in batches:
        threads_list_run(b)
    print("All chromosomes aligned")


# ============================================================================
# LOG PROCESSING FUNCTIONS (SERIAL)
# ============================================================================

def parse_conform_gt_log(log_file):
    """Parse Conform-GT log file into DataFrame."""
    df = pd.read_csv(log_file, sep='\t')
    return df


def save_full_csv(log_file):
    """Save full log as CSV."""
    df = parse_conform_gt_log(log_file)
    csv_file = log_file.replace('.log', '.csv')
    df.to_csv(csv_file, index=False)


def save_retain_from_study_csv(log_file):
    """Save filtered CSV (REMOVED variants to retain)."""
    df = parse_conform_gt_log(log_file)
    
    retain_info = ["NOT_IN_REFERENCE", "MULTIPLE_REF_MATCHES", "DUPLICATE_MARKER"]
    filtered = df[(df['SUMMARY'] == 'REMOVED') & (df['INFO'].isin(retain_info))]
    
    csv_file = log_file.replace('.log', '_RetainFromStudy.csv')
    filtered.to_csv(csv_file, index=False)
    
    return csv_file


def process_all_logs_serial(output_dir, output_prefix):
    """Process all log files serially."""
    
    print("Processing logs...")
    retained_csv_files = []
    
    for chromosome in range(1, 27):
        log_file = os.path.join(output_dir, f"{output_prefix}_{chromosome}_aligned.log")
        
        if not os.path.isfile(log_file):
            print(f"Chr {chromosome}: Log file not found")
            continue
        
        save_full_csv(log_file)
        retain_csv = save_retain_from_study_csv(log_file)
        retained_csv_files.append(retain_csv)
    
    print("Logs processed")
    return retained_csv_files


# ============================================================================
# CONSOLIDATION FUNCTIONS
# ============================================================================

def consolidate_retain_csvs(retained_csv_files, output_dir, output_prefix):
    """Merge all _RetainFromStudy CSVs into one."""
    
    if not retained_csv_files:
        print("No _RetainFromStudy CSV files found")
        return None
    
    # Read and concatenate
    dfs = [pd.read_csv(f) for f in retained_csv_files]
    consolidated_df = pd.concat(dfs, ignore_index=True)
    
    # Remove duplicates
    consolidated_df = consolidated_df.drop_duplicates(subset=['ID', 'POS'], keep='first')
    
    # Save
    output_csv = os.path.join(output_dir, f"{output_prefix}_AllSNPs_ToBe_RetainedFromStudy.csv")
    consolidated_df.to_csv(output_csv, index=False)
    print(f"Consolidated {len(retained_csv_files)} CSV files into master CSV")
    
    return output_csv


def create_extract_file(consolidated_csv, output_dir, output_prefix):
    """Create PLINK extract file from consolidated CSV."""
    
    if not os.path.isfile(consolidated_csv):
        print(f"Consolidated CSV not found")
        return
    
    df = pd.read_csv(consolidated_csv)
    snp_ids = df['ID'].dropna().unique()
    
    extract_file = os.path.join(output_dir, f"{output_prefix}_AllSNPs_ToBe_RetainedFromStudy.Extract-File.txt")
    
    with open(extract_file, 'w') as f:
        for snp_id in snp_ids:
            f.write(f"{snp_id}\n")
    
    print(f"Extract file created with {len(snp_ids)} SNPs")


# ============================================================================
# MAIN PROGRAM
# ============================================================================

def main():
    print("Conform-GT Alignment Pipeline\n")
    
    # Get user input
    plink_prefix = input("PLINK BED prefix: ").strip()
    reference_vcf = input("Reference VCF path: ").strip()
    conform_gt_jar = input("Conform-GT JAR path: ").strip()
    plinkPath = input("PLINK 1.9 path: ").strip()
    plinkPath2 = input("PLINK 2 path: ").strip()
    output_dir = input("Output directory: ").strip()
    output_prefix = input("Output prefix: ").strip()
    max_workers_input = input("Max parallel workers (default 5): ").strip()
    max_workers = int(max_workers_input) if max_workers_input else 5
    
    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Validate inputs
    if not os.path.isfile(f"{plink_prefix}.bed"):
        print("Error: PLINK .bed file not found")
        return
    if not os.path.isfile(reference_vcf):
        print("Error: Reference VCF not found")
        return
    if not os.path.isfile(conform_gt_jar):
        print("Error: Conform-GT JAR not found")
        return
    
    # Stage 1: Convert PLINK BED to VCF
    print("\nStage 1: Converting PLINK BED to VCF")
    target_vcf_path = os.path.join(output_dir, f"{output_prefix}_StudyInputFile_forAlignment")
    convert_plink_to_vcf(plink_prefix, target_vcf_path, plinkPath, plinkPath2)
    target_vcf_path += ".vcf.gz"
    # Stage 2: Run Conform-GT alignment (parallel batches)
    print("\nStage 2: Running Conform-GT alignment")
    run_alignment_parallel(reference_vcf, target_vcf_path, conform_gt_jar,
                          output_dir, output_prefix, max_workers)
    
    # Stage 3: Process logs (serial)
    print("\nStage 3: Processing logs")
    retained_csv_files = process_all_logs_serial(output_dir, output_prefix)
    
    # Stage 4: Consolidation
    print("\nStage 4: Consolidating results")
    consolidated_csv = consolidate_retain_csvs(retained_csv_files, output_dir, output_prefix)
    
    if consolidated_csv:
        create_extract_file(consolidated_csv, output_dir, output_prefix)
    
    # Done
    print("\nPipeline complete")
    print(f"Output: {output_dir}")


if __name__ == '__main__':
    main()
