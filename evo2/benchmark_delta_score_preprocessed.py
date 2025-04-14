#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import pandas as pd
import torch
from evo2 import Evo2
from tqdm import tqdm

coding_df = pd.read_csv("./hg38_dataset/Preprocessed_coding.csv")  # Contains coding variants
noncoding_df = pd.read_csv("./hg38_dataset/Preprocessed_noncoding.csv")  # Contains noncoding variants

# Load the evo2 model
model = Evo2('evo2_7b')
device = 'cuda:0' if torch.cuda.is_available() else 'cpu'

def calculate_delta_scores(df, output_path):
    """Calculate delta scores using paper's exact method"""
    # Store sequences and indexes
    ref_seqs = []
    ref_seq_to_index = {}
    ref_seq_indexes = []
    var_seqs = []
    
    # Parse sequences
    for _, row in tqdm(df.iterrows(), total=len(df), desc="Preparing sequences"):
        ref_seq = row['ref_sequence']
        var_seq = row['alt_sequence']
        
        # Track unique reference sequences
        if ref_seq not in ref_seq_to_index:
            ref_seq_to_index[ref_seq] = len(ref_seqs)
            ref_seqs.append(ref_seq)
        
        ref_seq_indexes.append(ref_seq_to_index[ref_seq])
        var_seqs.append(var_seq)
    
    # Convert to numpy arrays
    ref_seq_indexes = np.array(ref_seq_indexes)
    
    # Score sequences (batched automatically)
    print("Scoring reference sequences...")
    ref_scores = model.score_sequences(ref_seqs)
    
    print("Scoring variant sequences...")
    var_scores = model.score_sequences(var_seqs)
    
    # Calculate delta scores
    delta_scores = np.array(var_scores) - np.array(ref_scores)[ref_seq_indexes]
    
    # Add to dataframe and save
    df['delta_score'] = delta_scores
    df[['chrom', 'pos', 'delta_score']].to_csv(output_path, index=False)
    print(f"Saved delta scores to {output_path}")

# Process both datasets
calculate_delta_scores(coding_df, "./hg38_dataset/DeltaScores_coding.csv")
calculate_delta_scores(noncoding_df, "./hg38_dataset/DeltaScores_noncoding.csv")