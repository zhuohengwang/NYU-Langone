#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np

from Bio import SeqIO
from tqdm import tqdm


# In[2]:


benchmark_path = r'./benchmarks/ClinVar_gnomAD_benchmark_with_predictions.csv'
dataset_path = r'./uniprot_dataset/uniprot_sprot.fasta'
preprocessed_dataset_path = r'./uniprot_dataset/preprocessed_uniprot_sprot.fasta'
preprocessed_score_path = r'./uniprot_dataset/preprocessed_uniprot_sprot_esm_scores.csv'


# In[3]:


benchmark_df = pd.read_csv(benchmark_path)
benchmark_df = benchmark_df[['uniprot_id', 'aa_change', 'ESM1b_score']]

uniprot_id_list = list(dict.fromkeys(benchmark_df['uniprot_id'].tolist()))
benchmark_df.head()


# In[4]:


with open(preprocessed_dataset_path, "w") as out_file:
    records_dict = {record.id.split('|')[1]: record for record in SeqIO.parse(dataset_path, "fasta")}
    
    for uniprot_id in tqdm(uniprot_id_list):
        if uniprot_id in records_dict:
            record = records_dict[uniprot_id]
            record.id = uniprot_id
            record.description = uniprot_id
            SeqIO.write(record, out_file, "fasta")
    del records_dict

preprocessed_dataset_length = sum(1 for _ in SeqIO.parse(preprocessed_dataset_path, "fasta"))
print(preprocessed_dataset_length)


# In[5]:


#!python esm_score_missense_mutations.py --input-fasta-file ./uniprot_dataset/preprocessed_uniprot_sprot.fasta --output-csv-file ./uniprot_dataset/preprocessed_uniprot_sprot_esm_scores.csv


# In[6]:


preprocessed_dataset_df = pd.read_csv(preprocessed_score_path)

filtered_df = preprocessed_dataset_df[preprocessed_dataset_df.set_index(['seq_id', 'mut_name']).index.isin(benchmark_df.set_index(['uniprot_id', 'aa_change']).index)]

benchmark_df.to_csv(r'./uniprot_dataset/benchmark_df.csv')
filtered_df.to_csv(r'./uniprot_dataset/filtered_df.csv')


# In[ ]:




