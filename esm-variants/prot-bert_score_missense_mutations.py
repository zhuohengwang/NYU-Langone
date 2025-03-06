#!/usr/bin/env python
# coding: utf-8

# In[1]:


import argparse
import numpy as np
import pandas as pd

from transformers import BertModel, BertTokenizer
import torch
import torch.nn as nn
import torch.nn.functional as F

from tqdm import tqdm
from Bio import SeqIO
import time

# Define AAorder
AAorder = ['K','R','H','E','D','N','Q','T','S','C','G','A','V','L','I','M','P','Y','F','W']

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

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

# Get logits for the sequences
def get_logits(seq, model, tokenizer, device=0):
    inputs = tokenizer(seq, return_tensors="pt", padding=True, truncation=True).to(device)
    model.to(device)
    
    with torch.no_grad():
        logits = model(**inputs).last_hidden_state
        
    return F.softmax(logits, dim=-1).cpu().numpy()

# Compute LLR scores for a dataset of sequences
def get_wt_LLR(input_df, model, tokenizer, device=0, silent=False):
    genes = input_df.id.values
    LLRs = []
    input_df_ids = []
    
    for gname in tqdm(genes, disable=silent):
        seq_length = input_df[input_df.id == gname].length.values[0]
        
        if seq_length <= 1022:
            dt = [(gname+'_WT',input_df[input_df.id==gname].seq.values[0])]
            
            logits = get_logits(dt, model, tokenizer, device)
            WTlogits = pd.DataFrame(logits[0, 1:-1, :], columns=AAorder, index=list(input_df[input_df.id==gname].seq.values[0])).T.iloc[4:24].loc[AAorder]
            WTlogits.columns = [j + ' ' + str(i + 1) for i, j in enumerate(WTlogits.columns)]

            # skip unrecognized amino acid
            split_columns = [i.split(' ')[0] for i in WTlogits.columns]
            if not all(element in AAorder for element in split_columns):
                continue
            
            # Normalize by the diagonal (as described in the original code)
            wt_norm = np.diag(WTlogits.loc[[i.split(' ')[0] for i in WTlogits.columns]])
            LLR = WTlogits - wt_norm
            LLRs.append(LLR)
            input_df_ids.append(gname)
        else:
            # Handle sequences longer than 1022 characters (tiling)
            long_seq = input_df[input_df.id == gname].seq.values[0]
            ints, M, M_norm = get_intervals_and_weights(len(long_seq), min_overlap=512, max_len=1022, s=20)
            
            dt = [(gname + '_WT_' + str(i), ''.join(list(np.array(list(long_seq))[idx]))) for i, idx in enumerate(ints)]
            logit_parts = []
            for dt_ in chunks(dt, 20):
                batch_labels, batch_strs, batch_tokens = tokenizer(dt_, return_tensors="pt", padding=True, truncation=True).to(device)
                with torch.no_grad():
                    results_ = model(batch_tokens).last_hidden_state.cpu().numpy()
                for i in range(results_.shape[0]):
                    logit_parts.append(results_[i, :, :])

            logits_full = np.zeros((len(long_seq), 33))
            for i in range(len(ints)):
                logit = np.zeros((len(long_seq), 33))
                logit[ints[i]] = logit_parts[i].copy()
                logit = np.multiply(logit.T, M_norm[i, :]).T
                logits_full += logit

            WTlogits = pd.DataFrame(logits_full, columns=AAorder, index=list(long_seq)).T.iloc[4:24].loc[AAorder]
            WTlogits.columns = [j.split('.')[0] + ' ' + str(i + 1) for i, j in enumerate(WTlogits.columns)]
            
            # Skip unrecognized amino acid
            split_columns = [i.split(' ')[0] for i in WTlogits.columns]
            if not all(element in AAorder for element in split_columns):
                continue

            wt_norm = np.diag(WTlogits.loc[[i.split(' ')[0] for i in WTlogits.columns]])
            LLR = WTlogits - wt_norm
            LLRs.append(LLR)
            input_df_ids.append(gname)

    return input_df_ids, LLRs

def chop(L,min_overlap=511,max_len=1022):
    return L[max_len-min_overlap:-max_len+min_overlap]

def intervals(L,min_overlap=511,max_len=1022,parts=None):
    if parts is None: parts = []
    #print('L:',len(L))
    #print(len(parts))
    if len(L)<=max_len:
        if parts[-2][-1]-parts[-1][0]<min_overlap:
      #print('DIFF:',parts[-2][-1]-parts[-1][0])
            return parts+[np.arange(L[int(len(L)/2)]-int(max_len/2),L[int(len(L)/2)]+int(max_len/2)) ]
        else:
            return parts
    else:
        parts+=[L[:max_len],L[-max_len:]]
        L=chop(L,min_overlap,max_len)
        return intervals(L,min_overlap,max_len,parts=parts)

def get_intervals_and_weights(seq_len,min_overlap=511,max_len=1022,s=16):
    ints=intervals(np.arange(seq_len),min_overlap=min_overlap,max_len=max_len)
    ## sort intervals
    ints = [ints[i] for i in np.argsort([i[0] for i in ints])]
    
    a=int(np.round(min_overlap/2))
    t=np.arange(max_len)
    
    f=np.ones(max_len)
    f[:a] = 1 / (1 + np.exp(-(t[:a]-a/2)/s))
    f[max_len-a:] = 1 / (1 + np.exp((t[:a]-a/2)/s))
    
    f0=np.ones(max_len)
    f0[max_len-a:] = 1 / (1 + np.exp((t[:a]-a/2)/s))
    
    fn=np.ones(max_len)
    fn[:a] = 1 / (1 + np.exp(-(t[:a]-a/2)/s))
    
    filt=[f0]+[f for i in ints[1:-1]]+[fn]
    M = np.zeros((len(ints),seq_len))
    for k,i in enumerate(ints):
        M[k,i] = filt[k]
    M_norm = M/M.sum(0)
    return (ints, M, M_norm)

# Function to calculate PLL for a sequence
def get_PLL(seq, model, tokenizer, device=0):
    logits = get_logits(seq, model, tokenizer, device)
    idx = [tokenizer.convert_tokens_to_ids(i) for i in seq]
    return np.sum(np.diag(logits[0, :, idx]))

def main(args):
    """
    Execute the main script logic.
    
    Parameters:
    args (Namespace): Arguments parsed from command line input.
    """

    start_time = time.time()

    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    print('Using {}.'.format('GPU' if device == 'cuda' else 'CPU (this may be much slower)'))

    input_df = fasta_to_dataframe(args.input_fasta_file)

    print('Loading the model ({})...'.format(args.model_name))
    tokenizer = BertTokenizer.from_pretrained(args.model_name, do_lower_case=False)
    model = BertModel.from_pretrained(args.model_name)

    print('Invoking the model...')
    input_df_ids, LLRs = get_wt_LLR(input_df, model, tokenizer, device)

    print('Saving results...')
    results = []
    for seq_id, LLR in zip(input_df_ids, LLRs):
        raw_seq_results = LLR.transpose().stack().reset_index().rename(columns = {'level_0': 'wt_aa_and_pos', 'level_1': 'mut_aa', 0: 'esm_score'})
        seq_results = pd.DataFrame({'seq_id': seq_id, 'mut_name': raw_seq_results['wt_aa_and_pos'].str.replace(' ', '') + raw_seq_results['mut_aa'], 'esm_score': raw_seq_results['esm_score']})
        results.append(seq_results)
        
    results = pd.concat(results).reset_index(drop=True)
    results.to_csv(args.output_csv_file, index=False)

    print(f'Time used = {(time.time() - start_time):.2f}s')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compute ESM effect scores for all possible missense variants in a set of protein sequences.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('--input-fasta-file', dest='input_fasta_file', required=True, metavar='/path/to/input_aa_seqs.fasta', help='Path to the input FASTA file with the amino-acid sequences to calculate ESM scores for.')
    parser.add_argument('--output-csv-file', dest='output_csv_file', required=True, metavar='./esm_missense_effect_scores.csv', help='Path to save the output CSV file.')
    parser.add_argument('--model-name', dest='model_name', metavar='model_name', help='Name of the protBERT model to use.')

    args = parser.parse_args()

    main(args)


# In[ ]:




