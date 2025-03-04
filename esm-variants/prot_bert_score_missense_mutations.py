import argparse
import pandas as pd
from Bio import SeqIO
import time
from transformers import AutoModel, AutoTokenizer

import torch
from tqdm import tqdm

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

def get_wt_LLR(input_df, model, tokenizer, device=0, silent=False): 
    AAorder = ['K', 'R', 'H', 'E', 'D', 'N', 'Q', 'T', 'S', 'C', 'G', 'A', 'V', 'L', 'I', 'M', 'P', 'Y', 'F', 'W']
    
    genes = input_df.id.values
    LLRs = []
    input_df_ids = []

    for gname in tqdm(genes, disable=silent):
        seq_entry = input_df[input_df.id == gname]
        seq_length = seq_entry.length.values[0]
        seq = seq_entry.seq.values[0]

        # Preprocess sequence for ProtBERT with fixed length (20 tokens)
        formatted_seq = " ".join(seq)  # Add spaces between amino acids
        tokens = tokenizer(formatted_seq, return_tensors="pt", padding='max_length', truncation=True, max_length=20).to(device)

        with torch.no_grad():
            outputs = model(**tokens)
            embeddings = outputs.last_hidden_state  # Extract embeddings

        # Extract relevant logits, ensuring the embeddings match the sequence length
        embeddings = embeddings.squeeze(0)[1:-1, :]  # Remove [CLS] and [SEP] tokens

        # Adjust number of tokens: if length is less than 20, pad with dummy values
        num_tokens = embeddings.shape[0]
        if num_tokens < 20:
            padding = torch.zeros((20 - num_tokens, embeddings.shape[1]), device=embeddings.device)
            embeddings = torch.cat([embeddings, padding], dim=0)
        
        # Construct DataFrame using sequence as index (fixed length, max 20)
        WTlogits = pd.DataFrame(embeddings.cpu().numpy().T, index=AAorder[:num_tokens])  # Use actual sequence length
        WTlogits.columns = [f"{aa} {i+1}" for i, aa in enumerate(seq[:num_tokens])]  # Match the correct number of tokens

        split_columns = [col.split(' ')[0] for col in WTlogits.columns]
        if not all(aa in AAorder for aa in split_columns):
            continue

        # Compute log-likelihood ratio (LLR)
        wt_norm = np.diag(WTlogits.loc[WTlogits.index])
        LLR = WTlogits - wt_norm

        LLRs.append(LLR)
        input_df_ids.append(gname)

    return input_df_ids, LLRs

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
    tokenizer = AutoTokenizer.from_pretrained(args.model_name)
    model = AutoModel.from_pretrained(args.model_name).to(device)

    print('Invoking the model...')
    input_df_ids, LLRs = get_wt_LLR(input_df, model, tokenizer, device=device)

    print('Saving results...')
    results = []
    for seq_id, LLR in zip(input_df_ids, LLRs):
        raw_seq_results = (
            LLR.transpose()
            .stack()
            .reset_index()
            .rename(columns={'level_0': 'wt_aa_and_pos', 'level_1': 'mut_aa', 0: 'bert_score'})
        )
        seq_results = pd.DataFrame({
            'seq_id': seq_id,
            'mut_name': raw_seq_results['wt_aa_and_pos'].str.replace(' ', '') + raw_seq_results['mut_aa'],
            'bert_score': raw_seq_results['bert_score']
        })
        results.append(seq_results)

    final_results = pd.concat(results, ignore_index=True)
    final_results.to_csv("protbert_results.csv", index=False)

    print(f'Time used = {(time.time() - start_time):.2f}s')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compute ProtBERT embeddings for protein sequences.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('--input-fasta-file', dest='input_fasta_file', required=True, metavar='/path/to/input_aa_seqs.fasta', help='Path to the input FASTA file with the amino-acid sequences to process.')
    parser.add_argument('--output-csv-file', dest='output_csv_file', required=True, metavar='./protbert_embeddings.csv', help='Path to save the output CSV file.')
    parser.add_argument('--model-name', dest='model_name', default='Rostlab/prot_bert', metavar='MODEL_NAME', help='Name of the ProtBERT model to use. Default is Rostlab/prot_bert.')

    args = parser.parse_args()

    main(args)
