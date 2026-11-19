import os
import pandas as pd
import numpy as np

def process_vcf_files(input_dir, skiprows, output_combined_csv):
    """
    Process all VCF.gz files in a directory and extract specific columns.
    
    Parameters:
    input_dir (str): Directory containing .vcf.gz files
    skiprows (int): Number of lines to skip at the start of each VCF file
    output_combined_csv (str): Full path for the combined output CSV file
    """
    # Get list of all .vcf.gz files in the directory
    vcf_files = [f for f in os.listdir(input_dir) if f.endswith('.vcf.gz')]
    
    if not vcf_files:
        print("No .vcf.gz files found in the specified directory.")
        return
    
    csv_files = []
    
    for vcf_file in vcf_files:
        print(f"Processing {vcf_file}...")
        vcf_path = os.path.join(input_dir, vcf_file)
        
        try:
            # Read the VCF file fully with specified columns and skiprows
            df = pd.read_csv(vcf_path, sep='\t', skiprows=skiprows, compression='gzip', 
                           usecols=[0,1,2,3,4,5,6,7,8,9])
            
            # Extract DR2 and IMP from INFO column
            def extract_info(info_str):
                info_dict = {}
                if pd.isna(info_str):
                    return info_dict
                
                for item in str(info_str).split(';'):
                    if '=' in item:
                        key, value = item.split('=', 1)
                        info_dict[key] = value
                    else:
                        info_dict[item] = True
                return info_dict
            
            def parse_dr2_value(dr2_str):
                """
                Parse DR2 value which can be single value or comma-separated values.
                Returns the first value for multi-allelic sites, or None if invalid.
                """
                if not dr2_str or pd.isna(dr2_str):
                    return None
                
                try:
                    # Handle comma-separated values (multi-allelic)
                    if ',' in str(dr2_str):
                        # Take the first DR2 value for multi-allelic sites
                        first_value = str(dr2_str).split(',')[0].strip()
                        return float(first_value)
                    else:
                        # Single value
                        return float(dr2_str)
                except (ValueError, AttributeError):
                    return None
            
            # Apply extraction to INFO column
            info_data = df['INFO'].apply(extract_info)
            
            # Create DR2 column - extract and parse DR2 value
            df['DR2'] = info_data.apply(lambda x: parse_dr2_value(x.get('DR2')) if 'DR2' in x else None)
            
            # Create IMP_GNT column - 'IMP' if IMP flag is present, otherwise 'GT'
            df['IMP_GNT'] = info_data.apply(lambda x: 'IMP' if 'IMP' in x else 'GT')
            
            # Select required columns
            df_out = df[['#CHROM', 'POS', 'ID', 'QUAL', 'DR2', 'IMP_GNT']]
            
            # Save to CSV with _info.csv suffix
            csv_file = vcf_file.replace('.vcf.gz', '_info.csv')
            csv_path = os.path.join(input_dir, csv_file)
            df_out.to_csv(csv_path, index=False)
            csv_files.append(csv_path)
            print(f"Saved {csv_file}")
            
        except Exception as e:
            print(f"Error processing {vcf_file}: {str(e)}")
            continue
    
    # Combine all CSV files into one large CSV
    if csv_files:
        print("Combining all CSV files...")
        combined_df = pd.concat([pd.read_csv(f) for f in csv_files], ignore_index=True)
        combined_df.to_csv(output_combined_csv, index=False)
        print(f"Combined CSV saved as: {output_combined_csv}")
        print(f"Total rows in combined file: {len(combined_df)}")
    else:
        print("No CSV files were created to combine.")

def main():
    """
    Main function to get user inputs and run the processing.
    """
    print("VCF File Processor")
    print("=" * 50)
    
    # Get user inputs
    input_dir = input('Enter the directory path containing .vcf.gz files: ').strip()
    
    # Validate directory
    if not os.path.exists(input_dir):
        print("Error: Directory does not exist.")
        return
    
    skiprows = int(input('Enter the number of lines to skip at the start of each VCF file: '))
    output_combined_csv = input('Enter the full output path and filename for the combined CSV: ').strip()
    
    # Validate output directory
    output_dir = os.path.dirname(output_combined_csv)
    if output_dir and not os.path.exists(output_dir):
        print("Error: Output directory does not exist.")
        return
    
    # Process the files
    process_vcf_files(input_dir, skiprows, output_combined_csv)
    print('Processing complete.')

if __name__ == "__main__":
    main()
