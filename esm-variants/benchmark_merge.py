#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np


# In[2]:


benchmark_df_path = r'./uniprot_dataset/benchmark_df.csv'

reproduced_dict = {'esm1b_t33_650M_UR50S' : {'name' : 'ESM1b',
                                             'score_path' : r'./uniprot_dataset/ESM1b_reproduced_score.csv'
                                            },
                   'esm1v_t33_650M_UR90S_1' : {'name' : 'ESM1v-1',
                                               'score_path' : r'./uniprot_dataset/ESM1v-1_reproduced_score.csv'
                                              },
                   'esm2_t33_650M_UR50D' : {'name' : 'ESM2',
                                            'score_path' : r'./uniprot_dataset/ESM2_reproduced_score.csv'
                                           },
                   'prot-bert-Rostlab' : {'name' : 'ProtBERT',
                                          'score_path' : r'./uniprot_dataset/ProtBERT_reproduced_score.csv'
                                         }
                  }

reproduced_df_path = r'./uniprot_dataset/reproduced_df.csv'


# In[3]:


benchmark_df = pd.read_csv(benchmark_df_path)
benchmark_df.head()


# In[5]:


merge_df = None
for model_pretrained, model_dict in reproduced_dict.items():
    score_df = pd.read_csv(model_dict['score_path'])

    if 'esm_score' in score_df.columns:
        score_df = score_df.rename(columns={'esm_score': model_dict['name'] + '_score'})
    
    if 'ProtBERT_score' in score_df.columns:
        score_df = score_df.rename(columns={'ProtBERT_score': model_dict['name'] + '_score'})
    
    reproduced_df = score_df[
        score_df[['seq_id', 'mut_name']].apply(tuple, axis=1).isin(
            benchmark_df[['uniprot_id', 'aa_change']].apply(tuple, axis=1)
        )
    ]

    reproduced_df = reproduced_df.set_index(['seq_id', 'mut_name']).reindex(
        benchmark_df.set_index(['uniprot_id', 'aa_change']).index, 
        fill_value=None
    ).reset_index().rename(
        columns={'uniprot_id': 'seq_id', 'aa_change': 'mut_name'}
    )[['seq_id', 'mut_name', model_dict['name'] + '_score']].dropna(subset=[model_dict['name'] + '_score'])

    print(len(reproduced_df))
    if merge_df is None:
        merge_df = reproduced_df
    else:
        merge_df = pd.merge(merge_df, reproduced_df, on=['seq_id', 'mut_name'], how='inner')

# merge_df = reduce(lambda left, right: pd.merge(left, right, on=['seq_id', 'mut_name'], how='inner'), model_df_list)
merge_df.to_csv(reproduced_df_path, index=False)


# In[ ]:




