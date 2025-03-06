import argparse
import time
import numpy as np
import pandas as pd
from Bio import SeqIO
from proteinbert import load_pretrained_model
from proteinbert.conv_and_global_attention_model import get_model_with_hidden_layers_as_outputs
from tqdm import tqdm

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

# Function to compute global representations (proxy for log-likelihood)
def get_log_likelihood(seq, model, input_encoder):
    encoded_x = input_encoder.encode_X([seq], seq_len)  # Encode input
    _, global_representations = model.predict(encoded_x, batch_size=32, verbose=0)  # Get global representation
    
    # Apply softmax
    softmax_probs = tf.nn.softmax(global_representations, axis=-1)
    
    # Compute log-likelihood (log of probabilities)
    log_likelihood = np.log(softmax_probs + 1e-10)  # Small epsilon to avoid log(0)
    
    return np.mean(log_likelihood)

def main(args):
    """
    Execute the main script logic.
    
    Parameters:
    args (Namespace): Arguments parsed from command line input.
    """

    start_time = time.time()

    pretrained_model_generator, input_encoder = load_pretrained_model()
    seq_len = 1020 + 2  # Set sequence length (adjust as needed)
    model = get_model_with_hidden_layers_as_outputs(pretrained_model_generator.create_model(seq_len))

    input_df = fasta_to_dataframe(args.input_fasta_file)

    genes = input_df.id.values
    results = []
    
    for gname in tqdm(genes):
        seq_entry = input_df[input_df.id == gname]
        wild_seq = seq_entry.seq.values[0]
    
        # Compute wild-type log-likelihood
        ll_wt = get_log_likelihood(wild_seq, model, input_encoder)
    
        # Generate all single-point mutations and compute LLR
        for i, original_res in enumerate(wild_seq):
            for new_res in AAorder:
                mutated_seq = wild_seq[:i] + new_res + wild_seq[i+1:]
                ll_mut = get_log_likelihood(mutated_seq, model, input_encoder)
                llr = ll_mut - ll_wt  # Compute log-likelihood ratio
                
                mutation_label = f"{original_res}{i+1}{new_res}"
                results.append([gname, mutation_label, llr])
    
    # Convert results to DataFrame
    df = pd.DataFrame(results, columns=["seq_id", "mut_name", "score"])
    df.to_csv(args.output_csv_file, index=False)
    
    print("Mutation LLR scores saved")

    print(f'Time used = {(time.time() - start_time):.2f}s')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compute ProtBERT embeddings for protein sequences.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('--input-fasta-file', dest='input_fasta_file', required=True, metavar='/path/to/input_aa_seqs.fasta', help='Path to the input FASTA file with the amino-acid sequences to process.')
    parser.add_argument('--output-csv-file', dest='output_csv_file', required=True, metavar='./protbert_embeddings.csv', help='Path to save the output CSV file.')
    parser.add_argument('--model-name', dest='model_name', metavar='MODEL_NAME', help='Name of the ProtBERT model to use.')

    args = parser.parse_args()

    main(args)




