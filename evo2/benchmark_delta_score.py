#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import torch
import pandas as pd
from tqdm import tqdm
from evo2 import Evo2

# Load the evo2 model
evo2_model = Evo2('evo2_7b')

def compute_delta_score(ref_sequence, alt_sequence, ref, alt):
    """
    Compute the delta score for a given variant.
    """
    # Tokenize sequences
    ref_input_ids = torch.tensor(
        evo2_model.tokenizer.tokenize(ref_sequence),
        dtype=torch.int,
    ).unsqueeze(0).to(device)
    
    alt_input_ids = torch.tensor(
        evo2_model.tokenizer.tokenize(alt_sequence),
        dtype=torch.int,
    ).unsqueeze(0).to(device)
    
    # Get logits
    with torch.no_grad():
        ref_outputs, _ = evo2_model(ref_input_ids)
        alt_outputs, _ = evo2_model(alt_input_ids)
    
    ref_logits = ref_outputs[0]
    alt_logits = alt_outputs[0]
    
    # Get token indices for ref and alt alleles
    ref_token = evo2_model.tokenizer.tokenize(ref)[0]
    alt_token = evo2_model.tokenizer.tokenize(alt)[0]
    
    # Compute log probabilities
    ref_log_prob = torch.log_softmax(ref_logits, dim=-1)[0, len(ref_sequence) // 2, ref_token]
    alt_log_prob = torch.log_softmax(alt_logits, dim=-1)[0, len(alt_sequence) // 2, alt_token]
    
    # Compute delta score
    delta_score = (alt_log_prob - ref_log_prob).item()
    return delta_score

def compute_delta_scores(preprocessed_csv, output_path):
    """
    Compute delta scores for preprocessed data and save to a CSV file.
    """
    # Load preprocessed data
    preprocessed_df = pd.read_csv(preprocessed_csv)
    
    # Compute delta scores
    results = []
    
    for _, row in tqdm(preprocessed_df.iterrows(), desc="Computing delta scores", total=len(preprocessed_df)):
        chrom = row["chrom"]
        pos = row["pos"]
        ref = row["ref"]
        alt = row["alt"]
        ref_sequence = row["ref_sequence"]
        alt_sequence = row["alt_sequence"]
        
        # Compute delta score
        delta_score = compute_delta_score(ref_sequence, alt_sequence, ref, alt)
        
        results.append({
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "delta_score": delta_score
        })
    
    # Convert results to DataFrame
    results_df = pd.DataFrame(results)
    
    # Save results to a CSV file
    results_df.to_csv(output_path, index=False)
    print(f"Delta scores computed. Results saved to '{output_path}'.")

# Compute delta scores for coding variants
compute_delta_scores("./hg38_dataset/Preprocessed_coding.csv", "./hg38_dataset/DeltaScores_coding.csv")

# Compute delta scores for noncoding variants
compute_delta_scores("./hg38_dataset/Preprocessed_noncoding.csv", "./hg38_dataset/DeltaScores_noncoding.csv")

