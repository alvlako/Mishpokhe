import numpy as np
import pandas as pd

import itertools
import random
import subprocess


# This script simulate mapping results used in py_mishpokhe2.py for testing

n_enries = 100

df_index = list(range(0, n_enries))
df_comment = list(range(0, n_enries))
df_sort_cat = list(range(0, n_enries))


target_ids = list(range(0, n_enries))
query_ids = list(range(0, round(n_enries/10))) + 2*random.sample(range(0, round(n_enries/10)),  round(n_enries/10/2))+ list(range(0, round(n_enries/10)))+ 4*random.choices(range(0, round(n_enries/10)), k = round(n_enries/10))+ 3*list(range(0, round(n_enries/10))) 

df_strand =  round(n_enries/10)*[0] + 2*random.choices(range(0, 1),  k = round(n_enries/10/2))+ round(n_enries/10)*[1] + 4*random.choices(range(0, 1), k = round(n_enries/10))+ 3*round(n_enries/10)*[1]

df_coord1 = list()
df_coord2 = list()
start = 0
width = 10
for i in range(0, n_enries):
    coord1 = random.randrange(start, width)
    df_coord1.append(coord1)
    coord2 = coord1 + random.randrange(500, 1000)
    df_coord2.append(coord2)
    start = coord2 + 1
    width = coord2 + 1000
    i = i + 1


def return_mapped_res():
    mapped_results = pd.DataFrame(
    {'index': df_index,
     'ID': target_ids,
     'coord1': df_coord1,
     'coord2': df_coord2,
     'strand': df_strand,
     'comment': df_comment,
     'sort_cat': df_sort_cat,
     'query_ID': query_ids,
    })
    return (mapped_results)

# Testing

import sys
import os
sys.path.append(os.path.abspath('/Users/Sasha/Desktop'))

import py_mishpokhe2 as msh

mapped_res = return_mapped_res()

cluster_matches = msh.find_clusters()
print(cluster_matches)

print(x)

significant_cluster_df_enriched = msh.update_scores_for_cluster_matches(cluster_matches)
print(significant_cluster_df_enriched)
    
stat_lambda, stat_K = msh.calculate_karlin_stat(significant_cluster_df_enriched)
msh.calculate_e_value(stat_lambda, stat_K, significant_cluster_df_enriched)

msh.set_strand_flip_penalty(cluster_matches)
sign_clusters_df = msh.set_strand_flip_penalty(cluster_matches)

msh.extract_proteins_cluster_neighborhood(sign_clusters_df)
msh.update_query_profiles()
msh.add_new_proteins()
msh.make_new_profiles(sign_clusters_df)

msh.initialize_new_prot_score()
msh.combine_all_profiles()
