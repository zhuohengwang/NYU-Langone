#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np

from Bio import SeqIO
from tqdm import tqdm

import subprocess
from pathlib import Path


# In[2]:


benchmark_path = r'./benchmarks/ClinVar_gnomAD_benchmark_with_predictions.csv'
dataset_path = r'./uniprot_dataset/uniprot_sprot.fasta'
benchmark_df_path = r'./uniprot_dataset/benchmark_df.csv'
preprocessed_dataset_path = r'./uniprot_dataset/preprocessed_uniprot_sprot.fasta'

# note: all models below have 33 layers and 650M params
reproduced_dict = {'esm1b_t33_650M_UR50S' : {'name' : 'ESM1b',
                                             'score_path' : r'./uniprot_dataset/ESM1b_reproduced_score.csv',
                                             'calculator' : r'esm_score_missense_mutations.py'
                                            },
                   'esm1v_t33_650M_UR90S_1' : {'name' : 'ESM1v-1',
                                               'score_path' : r'./uniprot_dataset/ESM1v-1_reproduced_score.csv',
                                               'calculator' : r'esm_score_missense_mutations.py'
                                              },
                   'esm2_t33_650M_UR50D' : {'name' : 'ESM2',
                                            'score_path' : r'./uniprot_dataset/ESM2_reproduced_score.csv',
                                            'calculator' : r'esm_score_missense_mutations.py'
                                           }
                  }


# In[3]:


benchmark_df = pd.read_csv(benchmark_path)

# filter out row that has >1 scoring for them
other_score_columns = benchmark_df.columns[6:25]
benchmark_df = benchmark_df[benchmark_df[other_score_columns].notna().sum(axis=1) > 1]
benchmark_df = benchmark_df[['uniprot_id', 'aa_change', 'ESM1b_score', 'clinvar_label']]

benchmark_df.to_csv(benchmark_df_path, index=False)
print(len(benchmark_df))

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


# In[ ]:


for model_pretrained, model_dict in reproduced_dict.items():
    if Path(model_dict['score_path']).is_file():
        print(f"{model_dict['name']} already has score")
        continue
    subprocess.run([
        "python", model_dict['calculator'],
        "--input-fasta-file", preprocessed_dataset_path,
        "--output-csv-file", model_dict['score_path'],
        "--model-name", model_pretrained,
    ], check=True)


# In[ ]:




