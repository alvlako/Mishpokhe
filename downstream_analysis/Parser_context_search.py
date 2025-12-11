# Here is the script to organize the results of context search

from cmath import inf
import ctypes as ct
from ctypes import *
import numba as nb
import numpy as np
import pandas as pd
#import plotly.express as px
#import plotly.figure_factory as ff
#import plotly.io as pio

import argparse
import ast
import copy
import glob
import itertools
import logging
import os
import re
import subprocess
import sys
from sys import stdout
import time
#import packaging
#import networkx as nx
#import matplotlib.pyplot as plt

#pd.set_option('display.max_columns', None)


# Except for the parsing, I also add now the ability to compare to eggnog-annots of 2 columns with ' - ' separator

# Here it is parsing the mishpokhe-ncbi annots
# At the end of the loop I will also make it to the format of 2 columns: protein of interest, annots. Here I create this dataframe
mishpokhe_ncbi_annots_df = pd.DataFrame(columns=['prot_id','mish_annots'])
i = 0
for file in glob.glob("*.faa_context_res1iter_res_2_iter_sign_clusters_enrich_stat_filtered_clu_filter_ncbi_prods"):
    print(file)
    #cluster_annots = pd.read_csv(file, dtype=None, sep=',', header = None,converters={0: lambda x: x.split(", ")})
    try:
        cluster_annots = pd.read_csv(file, dtype=None, sep='\t', header = None)
    except pd.errors.EmptyDataError:
        print(f'{file} HAS NO ANNOTS')
        continue
    #,converters={0: lambda x: x.split(", ")})
    #print(cluster_annots)
    cluster_filter_file = file[:file.find('ncbi_prods')-1]
    #print(cluster_filter_file)
    cluster_filter = pd.read_csv(cluster_filter_file, dtype=None, sep='\t', converters={'query_prots':ast.literal_eval})
    #print('cluster_filter', cluster_filter)
    query_prots_col = cluster_filter['query_prots']

    # Now let's find what correponds to the protein of interest
    # First I will make the 2d numpy arrays
    cluster_annots_arr = np.array(pd.DataFrame(cluster_annots[0].str.split(',').values.tolist()).values)
    #query_prots_arr = np.array(cluster_filter['query_prots'].tolist(),dtype='str')
    query_prots_arr = np.array(pd.DataFrame(cluster_filter['query_prots'].values.tolist()).values)
    #print('cluster_annots_arr', cluster_annots_arr[0][0])
    #print('query_prots_arr', query_prots_arr[0][0])
    #print(np.where(cluster_annots_arr == 'integrase'))
    #print(np.where(query_prots_arr == 'NC_016765.1_33'))
    interest_prot = file[:file.find('.faa_context_res1iter_res_2_iter_sign_clusters_enrich_stat_filtered_clu_filter_ncbi_prods')]
    print('interest_prot', interest_prot)
    print('file', file)
    interest_prot_ind = np.where(query_prots_arr == interest_prot)
    #print(interest_prot_ind)
    #print('cluster_annots_arr', cluster_annots_arr)
    #print('query_prots_arr', query_prots_arr)
    print('shapes', np.shape(cluster_annots_arr), np.shape(query_prots_arr))
    # Check where the number of proteins in the arrays differ (probably because NCBI-defined ORFs might differ from mine, therefore more/less proteins)
    #print(cluster_annots_arr)
    interest_prot_annots = cluster_annots_arr[interest_prot_ind] 
    #print('interest_prot_annots', interest_prot_annots)
    #unique, counts = np.unique(interest_prot_annots, return_counts=True)
    #print('interest_prot annots', dict(zip(unique, counts)))
    #print(cluster_annots_arr)
    nonflat_cluster_annots_list = cluster_annots[0].str.split(',').values.tolist()
    #print('nonflat_cluster_annots_list', nonflat_cluster_annots_list)
    flat_cluster_annots_list = [x for xs in nonflat_cluster_annots_list for x in xs]
    #print('flat_cluster_annots_list', flat_cluster_annots_list)
    unique1, counts1 = np.unique(np.array(flat_cluster_annots_list), return_counts=True)
    # The argsort is to make the descending order
    neighbourhood_counts_dict = dict(zip(unique1[np.argsort(-counts1)], counts1[np.argsort(-counts1)]))
    print('neighbourhood annots', neighbourhood_counts_dict)
    mishpokhe_ncbi_annots_df.loc[i] = [interest_prot, neighbourhood_counts_dict]
    i = i + 1

# I need to sort the dataframe by prot_id
mishpokhe_ncbi_annots_df_sorted = mishpokhe_ncbi_annots_df.sort_values(by=['prot_id']).reset_index(drop=True)
print(mishpokhe_ncbi_annots_df_sorted)

# Now let's read the eggnog annots
pmg1_eggnog_df = pd.read_csv('pmg1_eggnog', dtype=None, sep=' - ', header = None)
#print(pmg1_eggnog_df)
pmg1_eggnog_df.columns = ['prot_id','egg_annots']
# I have to sort them just like mishpokhe annots, because the sorting is not completely numerical and it is complicated to fix
pmg1_eggnog_df_sorted = pmg1_eggnog_df.sort_values(by=['prot_id']).reset_index(drop=True)
print(pmg1_eggnog_df_sorted)

# The problem is also that mishpokhe might annotate not all proteins so now there are less mishpokhe annots than rows in eggnong dataframe (even though many of the are empty). Let's try to intersect dataframes to get rid of that
merged_annots = pd.merge(pmg1_eggnog_df_sorted, mishpokhe_ncbi_annots_df_sorted, how='outer', on='prot_id')
print(merged_annots)

# Okay, let's now split rows to words and compare sets of them
# here i would create a list for indices of rows in the dataframe where mishpokhe-ncbi and eggnog annotations had word intersection (for humans to check)
intersect_ind = []
for i in merged_annots.index:
    print(i)
    print(merged_annots['egg_annots'][i])
    if merged_annots['egg_annots'][i] == None:
        continue
    if pd.isna(merged_annots['mish_annots'][i]):
        continue
    egg_set = set(merged_annots['egg_annots'][i].split(' '))
    print('egg_set', egg_set)
    print(merged_annots['mish_annots'][i])
    mish_string = ' '.join(list(merged_annots['mish_annots'][i].keys()))
    mish_set = set(mish_string.split(' '))
    print('mish_set', mish_set)
    # here i want to remove non-meaningful annotations
    egg_set.difference_update({'protein', 'phage', 'bacterial', 'bacteria', 'bacteriophage', 'family', 'hypothetical', 'Hypothetical'})
    mish_set.difference_update({'protein', 'phage', 'bacterial', 'bacteria', 'bacteriophage', 'family', 'hypothetical', 'Hypothetical'})
    print('egg_set', egg_set)
    print('mish_set', mish_set)
    # Let's check if intersection is NOT empty
    if not egg_set.intersection(mish_set) == False:
        intersect_ind.append(i)

prelim_similar_annots = merged_annots.iloc[intersect_ind]
print('prelim_similar_annots', prelim_similar_annots)

prelim_similar_annots.to_csv('pmg1_prelim_similar_annots', sep = '\t')
