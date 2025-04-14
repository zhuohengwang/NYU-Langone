#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pysam
from Bio import SeqIO
import torch
from evo2 import Evo2
from tqdm import tqdm
import pandas as pd

# Load the reference genome
reference_genome = SeqIO.index("./hg38_dataset/hg38.fa", "fasta")

# Load the evo2 model
evo2_model = Evo2('evo2_7b')
device = 'cuda:0' if torch.cuda.is_available() else 'cpu'

def get_sequence(chrom, pos, ref, alt, window_size=512):
    """
    Extract a sequence window centered at the variant position.
    """
    start = pos - 1 - (window_size // 2)  # Convert to 0-based index
    end = start + window_size
    
    try:
        ref_sequence = str(reference_genome[chrom].seq[start:end])
        alt_sequence = ref_sequence[:window_size // 2] + alt + ref_sequence[window_size // 2 + len(ref):]
        return ref_sequence, alt_sequence
    except KeyError:
        print(f"Warning: Chromosome '{chrom}' not found in reference genome. Skipping variant.")
        return None, None
    except IndexError:
        print(f"Warning: Position {pos} is out of bounds for chromosome '{chrom}'. Skipping variant.")
        return None, None

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

def calculate_and_save_delta_scores(vcf_path, output_path):
    """
    Calculate delta scores for all variants in the VCF file and save to a CSV file.
    """
    # Open the VCF file
    vcf_reader = pysam.VariantFile(vcf_path)
    total_records = sum(1 for _ in vcf_reader)
    vcf_reader = pysam.VariantFile(vcf_path)
    
    delta_scores = []
    chroms = []
    pos_list = []
    
    for record in tqdm(vcf_reader, total = total_records, desc="Processing variants"):
        chrom = "chr" + record.chrom
        pos = record.pos
        ref = record.ref
        alts = record.alts
        
        # Skip if there are multiple alternate alleles
        if len(alts) != 1:
            continue
        
        alt = alts[0]
        
        # Get reference and alternate sequences
        ref_sequence, alt_sequence = get_sequence(chrom, pos, ref, alt)
        
        # Skip if sequences could not be extracted
        if ref_sequence is None or alt_sequence is None:
            continue
        
        # Compute delta score
        delta_score = compute_delta_score(ref_sequence, alt_sequence, ref, alt)
        delta_scores.append(delta_score)

        chroms.append(record.chrom)
        pos_list.append(record.pos)
    
    # Save delta scores and labels to a CSV file
    results_df = pd.DataFrame({
        "chrom": chroms,
        "pos": pos_list,
        "delta_score": delta_scores
    })
    results_df.to_csv(output_path, index=False)
    print(f"Delta scores saved to '{output_path}'.")

# Calculate and save delta scores for noncoding variants
calculate_and_save_delta_scores(
    "./benchmark/ClinVar_NonCoding_SNV_PB.filtered.vcf", 
    "./hg38_dataset/DeltaScores_noncoding.csv"
)

# Calculate and save delta scores for coding variants
calculate_and_save_delta_scores(
    "./benchmark/ClinVar_Coding_SNV_PB.vcf", 
    "./hg38_dataset/DeltaScores_coding.csv"
)

