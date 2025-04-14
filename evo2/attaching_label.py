#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import csv
import pandas as pd

coding_scores = pd.read_csv("./hg38_dataset/DeltaScores_coding.csv")  # chrom, pos, delta_score
noncoding_scores = pd.read_csv("./hg38_dataset/DeltaScores_noncoding.csv")

def extract_labels(vcf_path):
    """Extract labels from unannotated trailing field by parsing raw VCF lines"""
    labels = {}
    
    with open(vcf_path, 'r') as vcf_file:
        reader = csv.reader(vcf_file, delimiter='\t')
        for row in reader:
            if row[0].startswith('#'):
                continue  # Skip header lines
                
            chrom = row[0]
            pos = int(row[1])
            trailing_value = float(row[7])  # 8th column (0-indexed as 7)
            
            if trailing_value in (0.0, 1.0):
                labels[(chrom, pos)] = int(trailing_value)
    
    print(f"Extracted {len(labels)} labels from trailing field")
    return labels

# Get labels for coding/noncoding variants
coding_labels = extract_labels("./benchmark/ClinVar_Coding_SNV_PB.vcf")
noncoding_labels = extract_labels("./benchmark/ClinVar_NonCoding_SNV_PB.filtered.vcf")

def add_labels(df, label_dict):
    """Add labels to DataFrame using (chrom, pos) as keys."""
    df['label'] = df.apply(lambda row: label_dict.get((str(row['chrom']).replace('chr', ''), row['pos']), None), axis=1)
    return df

# Add labels to both DataFrames
coding_scores = add_labels(coding_scores, coding_labels)
noncoding_scores = add_labels(noncoding_scores, noncoding_labels)

coding_scores.to_csv("./hg38_dataset/DeltaScores_coding_with_labels.csv", index=False)
noncoding_scores.to_csv("./hg38_dataset/DeltaScores_noncoding_with_labels.csv", index=False)
