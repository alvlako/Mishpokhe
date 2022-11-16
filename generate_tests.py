import numpy as np
import pandas as pd

import itertools
import random
import string
import subprocess


# This script simulate mapping results used in py_mishpokhe2.py for testing

n_query_entries = 25
n_target_entries = 300

test_query_prot_file = open('test_query_prot.faa', 'w')

aa_list = 'A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,O,S,U,T,W,Y,V,B,Z,J'.split(',')
print(aa_list)
queries_list = []
j = 0
fasta_header = 0
for i in range(round(n_query_entries/6)):
    seq_len = random.choice(range(40, 500))
    print('seq', seq_len)
    represent_seq = ''.join(random.choices(aa_list, k = seq_len)) + '*'
    print(len(represent_seq))
    queries_list.append(represent_seq)
    i = i + 1
    j = j + 1
    # If you want to label members by cluster (rep), it can be done on this step
    for m in range(round(n_query_entries/(n_query_entries/6))-1):
        fasta_header = fasta_header + 1
        test_query_prot_file.write('>SimulatedTEST'+ str(fasta_header) + '\n')

        add_seq_len = random.choice(range(2, 50))
        member_seq = represent_seq[0:random.choice(range(seq_len-20))] + ''.join(random.choices(aa_list, k = add_seq_len)) + '*'
        queries_list.append(member_seq)
        test_query_prot_file.write(member_seq + '\n')
        m = m + 1
        j = j + 1
        print('j',j)
        print('n', n_query_entries)
while j < n_query_entries:
    fasta_header = fasta_header + 1
    test_query_prot_file.write('>SimulatedTEST'+ str(fasta_header) + '\n')
    
    add_seq_len = random.choice(range(2, 50))
    member_seq = represent_seq[0:random.choice(range(seq_len-20))] + ''.join(random.choices(aa_list, k = add_seq_len)) + '*'
    queries_list.append(member_seq)
    test_query_prot_file.write(member_seq + '\n')
    j = j + 1
    print('j',j)
    print('n', n_query_entries)
    
test_query_prot_file.close()
print(queries_list)
print(len(queries_list))


print(x)


# Testing

import sys
import os
sys.path.append(os.path.abspath('/Users/Sasha/Desktop'))

import py_mishpokhe2 as msh

mapped_res = return_mapped_res()

cluster_matches = msh.find_clusters()
print(cluster_matches)



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
