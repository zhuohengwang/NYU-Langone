#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as mcolors


# In[2]:


AAorder = ['K','R','H','E','D','N','Q','T','S','C','G','A','V','L','I','M','P','Y','F','W']


# In[3]:


def load_missense_data(missense_path, mut_col, score_col):
    # Load missense mutation data
    missense_df = pd.read_csv(missense_path)

    # Extract position and mutant amino acid
    missense_df["original_aa"] = missense_df[mut_col].str[0]  # First letter is original AA
    missense_df["position"] = missense_df[mut_col].str.extract(r'(\d+)').astype(int)
    missense_df["mut_aa"] = missense_df[mut_col].str[-1]  # Get last character (mutated AA)

    # Pivot to create a heatmap-friendly format
    missense_pivot = missense_df.groupby(["mut_aa", "position"])[score_col].sum().unstack()
    missense_pivot = missense_pivot.reindex(AAorder)

    unique_positions = sorted(missense_df["position"].unique())
    position_labels = [f"{missense_df[missense_df['position'] == pos]['original_aa'].values[0]} {pos}"
                       for pos in unique_positions]
    print(len(position_labels))
    
    return missense_pivot, position_labels


# In[6]:


viridis_r = plt.colormaps.get_cmap("viridis_r")
darker_viridis_r = mcolors.ListedColormap(np.vstack((viridis_r(np.linspace(0, 0.85, 128)) ** 1.45,
                                                     viridis_r(np.linspace(0.85, 1, 128)) ** 1.05
                                                    )))


# In[7]:


def plot_heatmap(missense_pivot, position_labels, output_fig_path, figsize):
    # Combine both heatmaps
    full_df = missense_pivot

    plt.figure(figsize=figsize)
    ax = sns.heatmap(full_df, cmap=darker_viridis_r, vmax=0, vmin=-26.5, center=0, xticklabels=position_labels, yticklabels=1)

    cbar = ax.collections[0].colorbar
    cbar.ax.set_title("LLR")

    plt.xticks(rotation='vertical')
    plt.yticks(rotation='horizontal')

    plt.xlabel("Protein sequence")
    plt.ylabel("Amino acid change")
    plt.title("Missense Effects")

    plt.savefig(output_fig_path)
    plt.show()

'''
mut_col = "variant"
score_col = "score"
curr_gene = r'SPAST (SPAST) _ Q9UBP0'
missense_path = rf'./examples/{curr_gene}.csv'
output_fig_path = rf'./examples/{curr_gene}_gen.jpg'
missense_pivot, position_labels = load_missense_data(missense_path, mut_col, score_col)
plot_heatmap(missense_pivot, position_labels, output_fig_path, figsize = (120, 5))
'''


