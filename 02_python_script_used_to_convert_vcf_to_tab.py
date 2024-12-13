#!/usr/bin/env python3
"""
Author: g. ozan bozdag

This script processes VCF (Variant Call Format) files and extracts chromosome and position information.
It creates simplified TAB-delimited files containing only the CHROM and POS columns from the input VCF files.

Input: VCF files in the current directory
Output: TAB files with .tab extension containing CHROM and POS columns

Usage: Place this script in the directory containing VCF files and run it.
The script will process all .vcf files in the directory and create corresponding .tab files.
"""

import os

def extract_vcf_info(input_vcf, output_tab):
    """
    Extract chromosome and position information from VCF file and save to TAB format.
    
    Args:
        input_vcf (str): Path to input VCF file
        output_tab (str): Path to output TAB file
    """
    with open(input_vcf, 'r') as infile, open(output_tab, 'w') as outfile:
        for line in infile:
            # Skip VCF header lines starting with ##
            if line.startswith('##'):
                continue
            # Process column header line starting with #CHROM
            if line.startswith('#CHROM'):
                outfile.write("\t".join(line.strip().split('\t')[:2]) + '\n')
            # Process variant data lines
            else:
                parts = line.strip().split('\t')
                outfile.write('\t'.join(parts[:2]) + '\n')

if __name__ == "__main__":
    directory_path = "."  # Current directory. Modify as needed.
    # Process all VCF files in the directory
    for filename in os.listdir(directory_path):
        if filename.endswith('.vcf'):
            input_vcf = os.path.join(directory_path, filename)
            output_tab = os.path.join(directory_path, filename.replace('.vcf', '.tab'))
            extract_vcf_info(input_vcf, output_tab)
            print(f"Processed {input_vcf} -> {output_tab}")