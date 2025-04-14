import pandas as pd
import csv

def convert_csv_to_vcf(csv_path, vcf_output_path):
    """Convert CSV with tab-delimited VCF data to proper VCF file"""
    # Read the CSV file
    df = pd.read_csv(csv_path, sep='\t', quoting=csv.QUOTE_NONNUMERIC).rename(columns = {"#CHROM" : "CHROM"})
    df['POS'] = df['POS'].astype(int)
    df['ID'] = df['ID'].astype(int)
    print(df.head())
    valid_df = df[df['POS'] != -1].copy()
    filtered_count = len(df) - len(valid_df)

    binary_columns = [
        'coding_sequence', 'start_codon', 'stop_codon', 
        'five_prime_UTR', 'three_prime_UTR', 'intron',
        'promoter', 'splice', 'RNA_gene', 'other'
    ]
    
    # Create VCF header
    vcf_header = """##fileformat=VCFv4.2
##source=ClinVar
##INFO=<ID=ClinVar_annotation,Number=1,Type=String,Description="1: Pathogenic; 0:Benign">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
""" 
# + "\t".join(binary_columns) + "\n"

    # Write to VCF file
    with open(vcf_output_path, 'w') as vcf_file:
        vcf_file.write(vcf_header)
        
        for _, row in valid_df.iterrows():
#            binary_values = [str(row.get(col, 0)) for col in binary_columns]
#            tab_separated = '\t'.join(binary_values)

            # Reconstruct the VCF line
            vcf_line = (f"{row['CHROM']}\t{row['POS']}\t{row['ID']}\t"
                       f"{row['REF']}\t{row['ALT']}\t.\t.\t{row['INFO']}\n"
                       )
            vcf_file.write(vcf_line)
    
    print(f"Converted {len(valid_df)} variants (filtered {filtered_count} with POS=-1) to {vcf_output_path}")

# Example usage:
# convert_csv_to_vcf("./benchmark/ClinVar_NonCoding_splice_PB.csv", "./benchmark/ClinVar_NonCoding_splice_PB.vcf")
# convert_csv_to_vcf("./benchmark/ClinVar_NonCoding_NoncodingGene_PB.csv", "./benchmark/ClinVar_NonCoding_NoncodingGene_PB.vcf")

convert_csv_to_vcf("./benchmark/ClinVar_benchmark_all.csv", "./benchmark/ClinVar_benchmark_all.vcf")