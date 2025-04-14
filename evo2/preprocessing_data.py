#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pysam
from Bio import SeqIO
import pandas as pd
from tqdm import tqdm

# Load the reference genome
reference_genome = SeqIO.index("./hg38_dataset/hg38.fa", "fasta")

def get_sequence(chrom, pos, ref, alt, window_size=8192): #previously was 512
    """
    Extract a sequence window centered at the variant position.
    """
    start = pos - 1 - (window_size // 2)  # Convert to 0-based index
    end = start + window_size
    ref_sequence = str(reference_genome[chrom].seq[start:end])
    alt_sequence = ref_sequence[:window_size // 2] + alt + ref_sequence[window_size // 2 + len(ref):]
    return ref_sequence, alt_sequence

def preprocessing_data(vcf_reader, output_path, total_records):
    # Preprocess and save data
    preprocessed_data = []
    
    for record in tqdm(vcf_reader, total = total_records, desc="Processing variants"):
        chrom = f"chr{record.chrom}" if not record.chrom.startswith("chr") else record.chrom
        pos = record.pos
        ref = record.ref
        alt = record.alts[0]  # Assuming only one alternate allele
        
        # Get reference and alternate sequences
        ref_sequence, alt_sequence = get_sequence(chrom, pos, ref, alt)
        
        # Save metadata and sequences
        preprocessed_data.append({
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "ref_sequence": ref_sequence,
            "alt_sequence": alt_sequence
        })
    
    # Convert to DataFrame
    preprocessed_df = pd.DataFrame(preprocessed_data)
    
    # Save preprocessed data to a CSV file
    preprocessed_df.to_csv(output_path, index=False)

    print(f"Preprocessed {total_records} records. Data saved to {output_path}.")

#vcf_file = "./benchmark/ClinVar_Coding_SNV_PB.vcf"
#vcf_reader = pysam.VariantFile(vcf_file)
#total_records = sum(1 for _ in vcf_reader)
#vcf_reader = pysam.VariantFile(vcf_file)
#preprocessing_data(vcf_reader, "./hg38_dataset/Preprocessed_coding.csv", total_records)

#vcf_file = "./benchmark/ClinVar_NonCoding_SNV_PB.filtered.vcf"
#vcf_reader = pysam.VariantFile(vcf_file)
#total_records = sum(1 for _ in vcf_reader)
#vcf_reader = pysam.VariantFile(vcf_file)
#preprocessing_data(vcf_reader, "./hg38_dataset/Preprocessed_noncoding.csv", total_records)

vcf_file = "./benchmark/ClinVar_benchmark_all.vcf"
vcf_reader = pysam.VariantFile(vcf_file)
total_records = sum(1 for _ in vcf_reader)
vcf_reader = pysam.VariantFile(vcf_file)
preprocessing_data(vcf_reader, "./hg38_dataset/Preprocessed_all.csv", total_records)

