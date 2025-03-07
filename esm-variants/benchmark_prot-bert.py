#!/usr/bin/env python
# coding: utf-8

# In[1]:


import torch
import pandas as pd
from transformers import AutoModelForMaskedLM, AutoTokenizer
from Bio import SeqIO
from tqdm import tqdm  # For progress bar

# Load ProtBERT

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

model = AutoModelForMaskedLM.from_pretrained("Rostlab/prot_bert").to(device)
tokenizer = AutoTokenizer.from_pretrained("Rostlab/prot_bert")

# Convert FASTA file to DataFrame
def fasta_to_dataframe(fasta_file):
    """
    Convert a FASTA file into a Pandas DataFrame.
    
    Parameters:
    fasta_file (str): Path of the FASTA file.
    
    Returns:
    df (DataFrame): Pandas DataFrame containing id, gene, seq, and length columns.
    """
    records = list(SeqIO.parse(fasta_file, "fasta"))
    data = {
        "id": [record.id for record in records],
        "gene": [record.description.split()[1] if len(record.description.split()) > 1 else record.id for record in records],
        "seq": [str(record.seq) for record in records],
        "length": [len(record.seq) for record in records]
    }
    df = pd.DataFrame(data)
    return df

# Parse mutation string (e.g., "N660S")
def parse_mutation(mut_name):
    """Extract wild-type, position, and variant from mutation string."""
    wt_aa = mut_name[0]  # Wild-type amino acid (e.g., "N")
    position = int(mut_name[1:-1]) - 1  # Convert to 0-based index (e.g., 659)
    var_aa = mut_name[-1]  # Variant amino acid (e.g., "S")
    return wt_aa, position, var_aa

# Compute ProtBERT score for a mutation
def compute_protbert_score(sequence, seq_id, wt_aa, position, var_aa):
    """Compute LLR (log-likelihood ratio) for a mutation."""
    
    # Add spaces between amino acids for ProtBERT tokenizer
    spaced_sequence = " ".join(list(sequence))
    
    # Tokenize the sequence
    tokens = tokenizer(spaced_sequence, return_tensors="pt").to(device)
    input_ids = tokens["input_ids"]
    
    # Mask the position (account for [CLS] offset)
    mask_index = position + 1  # Skip [CLS]
    
    tokens["input_ids"][0, mask_index] = tokenizer.mask_token_id
    
    # Get logits and compute log probabilities
    with torch.no_grad():
        logits = model(**tokens).logits[0, mask_index]
    log_probs = torch.log_softmax(logits, dim=-1)
    
    # Calculate LLR
    aa_vocab = tokenizer.get_vocab()
    llr = log_probs[aa_vocab[var_aa]] - log_probs[aa_vocab[wt_aa]]
    return llr.item()

# Main function
def process_mutations(fasta_file, csv_file, output_file):
    """Process mutations and compute ProtBERT scores."""
    # Load sequences and mutations
    fasta_df = fasta_to_dataframe(fasta_file)
    mutations_df = pd.read_csv(csv_file)
    
    # Prepare output data
    results = []
    
    # Process each mutation with a progress bar
    for _, row in tqdm(mutations_df.iterrows(), total=len(mutations_df), desc="Processing mutations"):
        seq_id = row["uniprot_id"]
        mut_name = row["aa_change"]
        
        # Find the sequence in the FASTA DataFrame
        sequence_row = fasta_df[fasta_df["id"] == seq_id]
        if sequence_row.empty:
            print(f"Sequence {seq_id} not found in FASTA file. Skipping.")
            continue
        sequence = sequence_row["seq"].values[0]
        
        # Parse mutation
        wt_aa, position, var_aa = parse_mutation(mut_name)

        if position >= len(sequence) or sequence[position] != wt_aa:
            continue
            
        # Compute ProtBERT score
        score = compute_protbert_score(sequence, seq_id, wt_aa, position, var_aa)
        results.append({"seq_id": seq_id, "mut_name": mut_name, "ProtBERT_score": score})
    
    # Save results to CSV
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")

# Run the script
if __name__ == "__main__":
    fasta_file = r'./uniprot_dataset/preprocessed_uniprot_sprot.fasta'
    csv_file = r'./uniprot_dataset/benchmark_df.csv'
    output_file = r'./uniprot_dataset/ProtBERT_reproduced_score.csv'  # Use a different output file
    process_mutations(fasta_file, csv_file, output_file)


# In[ ]:




