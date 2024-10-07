from cmath import inf
import ctypes as ct
from ctypes import *
import numba as nb
import numpy as np
import pandas as pd

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

# version 2 accordingly to the written in the proposal

def arg_parser():
    if len(sys.argv) == 1:
        sys.exit("No arguments provided. You need to provide path to the files (see py_mishpokhe2.py -h)")
    parser = argparse.ArgumentParser(description="Mishpokhe: self-supervised discovery of functional clusters")
    parser.add_argument("-q", "--querydb",
     help="Provide path to query MMseqs2 database (protein, unordered)")
    parser.add_argument("-t", "--targetdb",
     help="Provide path to target MMseqs2 database (protein, ordered)")
    parser.add_argument("-r", "--res",
     help="Specify the name to be given to the results of the search")
    parser.add_argument("-tf", "--targetfa",
     help="Provide path to target sequence (fasta (linearized!))")
    parser.add_argument("-qf", "--queryfa",
     help="Provide path to query sequence (fasta)")
    parser.add_argument("-i", "--iter",
     help="Give the number of iterations, default is 1", default = 1)
    parser.add_argument("-s", "--singleton",
     help="Set to 1 if you have singleton queries, default is 0",  default=0)
    parser.add_argument("-u", "--evalfilteruse",
     help="Set to 1 if you want to use mishpokhe clusters e-value filter, default is 1",  default=1)
    parser.add_argument("-e", "--eval",
     help="Specify the e-value threshold, default is 100",  default=100)
    parser.add_argument("-f", "--frac_occ_min",
     help="Specify the threshold for the fraction of the cluster matches in which each sequence cluster occurs, default is 0",  default=0)
    parser.add_argument("-c", "--search_cov",
     help="Specify the coverage threshold for the mmseqs search",  default='0.8')
    parser.add_argument("-af", "--arc_filter",
     help="Set to 0 if you do NOT want to use architecture clustering and filtering, default is 1",  default=1)
    # make separate func?
    # CHECK is it ok to make global?
    global args
    args, unknown = parser.parse_known_args()
    print(vars(args))
    for argument in vars(args):
        arg_path = getattr(args, argument)
        if not os.path.exists((arg_path)):
            if argument not in ['res', 'iter', 'singleton', 'evalfilteruse', 'eval', 'frac_occ_min', 'search_cov', 'arc_filter']:
                sys.exit(f"{arg_path} not found")

class FilePath:
     """
     Reading paths to the files

     """
     def __init__(self, query_db, target_db, res):
        self.query_db = query_db
        self.target_db = target_db
        self.res = res


     @classmethod
     def get_path(self):
         while 1:
            try:
                query_db = args.querydb
                target_db = args.targetdb
                res = args.res
                return self(query_db,target_db, res)
            except:
                print('Invalid input!')
                continue


# add parameters adjustment option? error raising?
# deleting previous tmp and db files
def make_profiles():
    print('making query mmseqs profiles')
    # making query profile db from scratch is only for the first iteration
    subprocess.call(['mmseqs', 'cluster', files.query_db, files.query_db + '_clu', 'tmp', '-c', '0.8', '-e', '0.001'])
    
    subprocess.call(['mmseqs', 'result2msa', 
    files.query_db,
     files.query_db,
     files.query_db + '_clu',
     files.query_db + '_clu_msa',
     '--msa-format-mode', '4'])

    subprocess.call(['mmseqs', 'convertmsa', 
     files.query_db + '_clu_msa',
     files.query_db + '_clu_msa_db'])
    
    subprocess.call(['mmseqs', 'msa2profile', files.query_db + '_clu_msa_db', files.query_db + '_clu_msa_db_profile'])


# add options?
# CHECK if its right to search query profiles against target
# ADD besthit
def run_search():
    print('running mmseqs profile search')
    # extracting best hit per match with --max-accept 1
    # ASK Johannes if it is ok to use this without sensitivity iters (manual)
    # Ask Johannes if it is correct to search TARGET against Query
    # FIX, --max-accept 1 does not work for now
    # now the coverage and the e-value threshold is just as for the clustering
    search_cov = str(args.search_cov)
    subprocess.call(['mmseqs', 'search', 
    #files.query_db + '_clu' + '_rep' + '_profile',
    files.query_db + '_clu_msa_db_profile',
     files.target_db,
     files.res + '_prof_search',
     'tmp', '-a', '--mask', '0', '--comp-bias-corr', '0', '--max-seqs', '10000', '-c', search_cov, '-e', '0.001'])
    print(f'command on the {iter_counter} is mmseqs search  {files.query_db}_clu_msa_db_profile {files.target_db} {files.res}_prof_search tmp -a --mask 0 --comp-bias-corr 0 --max-seqs 10000 -c 0.8 -e 0.001')
    #, '--min-seq-id', '0.5'
    #subprocess.call(['mmseqs', 'search', files.target_db,
    #files.query_db + '_clu' + '_rep' + '_profile',
    # files.res + '_prof_search',
    # 'tmp', '-a'])
    #, '--max-accept', '1'])
    # convert to convenient format
    subprocess.call(['mmseqs', 'convertalis',
    # files.query_db + '_clu' + '_rep' + '_profile',
    files.query_db + '_clu_msa_db_profile',
     files.target_db, files.res + '_prof_search',
      files.res + '_prof_search' +'.m8'])
    #subprocess.call(['mmseqs', 'convertalis', 
    # files.target_db, files.query_db + '_clu' + '_rep' + '_profile', files.res + '_prof_search',
    #  files.res + '_prof_search' +'.m8'])
    # As i change the search direction, I have to swap the columns in the file to use the same 
    # code for other parts as before
    # THINK if should swap for non-coverted result file
    #with open("tmp_res", "w") as outfile:
    #    subprocess.call(["awk", "{ t = $1; $1 = $2; $2 = t; print; }", "OFS=\t", files.res + "_prof_search" +".m8"], stdout=outfile)
    #with open(files.res + "_prof_search" +".m8", "w") as outfile2:
    #    subprocess.call(["cat", "tmp_res"], stdout=outfile2)
    


# The n of prots in db matches the n of gff CDS CDS!!! entries
# FIX best query hit is needed
# CHECK if the mapping correct 
# CHECK if it quaranteed that indices in results file are same order as in lookup (also in order) 
# AND in _h and the order is always the same
# FIX !!!!!! some real ids are read as NaN
# currently it is search used the normal db for target
# CHECK if it is target proteins, not query
pd.set_option('display.max_columns', None)
class ResultsMapping:
     """
     Mapping results to their coordinates and strands. Produces file
     with 0 column with NaNs (FIX!!!), 1 column with coordinates and strand,
     ind column with the internal id of the protein

     """
     # the results file is currently m8 file from normal db search
     def __init__(self, search_result_file, target_db_lookup, target_db_h, res_map_to_header, ind_list):
        self.search_result_file = search_result_file
        self.target_db_lookup = target_db_lookup
        self.target_db_h = target_db_h
        self.res_map_to_header = res_map_to_header
        self.ind_list = ind_list


     @classmethod
     def map_target_to_coord(self, query_specific_thresholds):
        print('mapping results')
        #search_result_file = pd.read_csv(str(files.res + '_prof_search'), dtype={'int':'float'}, sep='\t')
        #search_result_file = pd.read_csv('/Users/Sasha/Documents/GitHub/mishpokhe_test/py2_multihit_res', dtype={'int':'float'}, sep='\t')
        #search_result_file = pd.read_csv('/Users/Sasha/Documents/GitHub/mishpokhe_test/py_norm_res_prof_search.m8', dtype={'str':'float'}, sep='\t', header = None)
        search_result_file = pd.read_csv(files.res + '_prof_search' +'.m8', dtype={'str':'float'}, sep='\t', header = None)
        logging.debug(f'search_result_file: {search_result_file}')
        #target_db_lookup = np.genfromtxt(str(files.target_db)+str(".lookup"), dtype = None, delimiter="\t", encoding=None)
        target_db_lookup = pd.read_csv(str(files.target_db)+str(".lookup"), dtype=None, sep='\t', header = None)
        target_db_h = pd.read_csv(str(files.target_db)+str("_h"), sep='\s+#\s+', header=None, engine='python')
        #print(search_result_file.iloc[:, 0])
        target_db_lookup.columns = ["ind_id", "ID", "file_id"]
        
        # RENAME ID to target_ID
        target_db_h.columns = ["ID", "coord1", "coord2", "strand","comment"]
        # CHECK if this pattern enough to remove?
        target_db_h["ID"] = target_db_h["ID"].str.replace('\x00', '')

        ##ind_list = search_result_file.iloc[:, 0].dropna().astype(int)
        #print(ind_list)
        # get target proteins real ids
        # FIX "unique"? is it ok to use?
        # assigning the best-matching query to a target (with min e-val)
        search_result_file.columns = ["query_id", "target_id","seq_ident", "score",
         "smth1", "smth2", "smth3", "smth4", "smth5", "smth6", "eval","score2"]
        #add_row = ['WIN93251.1','MW067000.1_7','0.981','207','8','0','1','207','1','207', '2.364000e-175', '501']
        #search_result_file.loc[len(search_result_file)] = add_row
        #add_row = ['SEWIN93251.1','MT006214.1_2','0.991','207','8','0','1','207','1','207', '2.364000e-178', '501']
        #search_result_file.loc[len(search_result_file)] = add_row
        search_result_file["eval"] = pd.to_numeric(search_result_file["eval"])
        logging.debug(f'search_result_file with cols: {search_result_file}')
        #print(search_result_file.groupby('target_id')["eval"].transform('min'))
        # that is to only get best-matching query
        search_result_file = search_result_file.loc[search_result_file.groupby('target_id')['eval'].idxmin()].reset_index(drop=True)

        # REMOVE unique? it's unique from the prev step
        real_id_list = pd.Series(pd.unique(search_result_file.iloc[:, 1]))
        logging.debug(f'real ids list: {real_id_list}')
        # map by 1st (0) column with real ids from search res
        # print(target_db_h.loc[target_db_h.iloc[:, 0].astype(str) == 'MT006214.1_1'])
        # dropping duplicates for these weird cases when prodigal proteins have duplicates with different comments
        tmp_res_map_to_header = target_db_h.loc[target_db_h.iloc[:, 0].astype(str).isin(real_id_list)].drop_duplicates(subset=['ID'])

        # NOTe - isin makes sorting, so the next lines to get to the original lines order
        # that is to link proper query ids to the other info
        tmp_res_map_to_header['sort_cat'] = pd.Categorical(tmp_res_map_to_header['ID'], categories=real_id_list, ordered=True)
        tmp_res_map_to_header.sort_values('sort_cat', inplace=True)
        tmp_res_map_to_header.reset_index(inplace=True)
        tmp_res_map_to_header['query_ID'] = search_result_file.iloc[:, 0]
        # now i need bit score for query specific match thresholds
        tmp_res_map_to_header['bit_score'] = search_result_file.iloc[:, 11]

        # sort it back again to iterate via by order in the next step
        # these ids are strings, need to be numeric for correct sorting
        
        def atoi(text):
            return int(text) if text.isdigit() else text

        def natural_keys(text):
            return [ atoi(c) for c in re.split(r'(\d+)', text) ]
        
        print(' ')
        print(' ')
        print(' ')
        print('start sorting')
        print(' ')

        # set is added to have only unique ids here in case there are duplicates in
        # target db
        human_ids = list(set(tmp_res_map_to_header['ID'].tolist()))
        human_ids.sort(key=natural_keys)
        #tmp_res_map_to_header['ID'] = human_ids.astype('category')
        tmp_res_map_to_header['ID_cat'] = pd.Categorical(tmp_res_map_to_header['ID'], categories=human_ids, ordered=True)
        #tmp_res_map_to_header.sort_values('ID_cat')
        
        tmp_res_map_to_header.sort_values('ID_cat', inplace=True)
        tmp_res_map_to_header.reset_index(inplace=True, drop=True)
        res_map_to_header = tmp_res_map_to_header.copy()
        #res_map_to_header = tmp_res_map_to_header.sort_values(by=['ID'])

        ##res_map_to_header['ind'] = ind_list.values
        
        # here i have to filter by match score query-specific thresholds
        print('res_map_to_header', res_map_to_header)
        print('query_specific_thresholds', query_specific_thresholds)
        res_map_to_header['match_thresholds'] = res_map_to_header['query_ID'].map(query_specific_thresholds)
        print('res_map_to_header2', res_map_to_header)
        res_map_to_header = res_map_to_header.drop(res_map_to_header[res_map_to_header['bit_score'] < res_map_to_header['match_thresholds']].index)
        print('res_map_to_header3', res_map_to_header)

        # sorting target db lookup to iterate correctly in find_clusters()
        # DOUBLE check?
        human_ids2 = list(set(target_db_lookup['ID'].tolist()))
        human_ids2.sort(key=natural_keys)
        target_db_lookup['ID_cat'] = pd.Categorical(target_db_lookup['ID'], categories=human_ids2, ordered=True)
        target_db_lookup.sort_values('ID_cat', inplace=True)
        target_db_lookup.reset_index(inplace=True, drop=True)
        target_db_lookup = target_db_lookup.drop(['ID_cat'], axis=1)

        # somehow target_db_h doesnt get always properly sorted in find_clusters, have to do it here
        
        target_db_h['ID_cat'] = pd.Categorical(target_db_h['ID'], categories=human_ids2, ordered=True)
        target_db_h.sort_values('ID_cat', inplace=True)
        target_db_h.reset_index(inplace=True, drop=True)
        target_db_h = target_db_h.drop(['ID_cat'], axis=1)


        # get target proteins real ids
        # CLEAN
        ind_list = real_id_list
        res_map_to_header.to_csv('res_map_to_header', sep = '\t')
        #print(res_map_to_header)
        return self(search_result_file, target_db_lookup, target_db_h, res_map_to_header, ind_list)
        #pass


# ?? is it correct that L is set not of all the target proteins
# but only those which were matched and got to results???
# FIX !!!!!! some real ids are read as NaN
# FIX ! to process not the initial results but the best hit results
# CHECK whether you need to change the scores
# DELETE header from pandas dataframes from ResultMapping class object
# Set the scores not to e-values column but to 0
# Now that is for the processing of norm search res, so indices from
# the mapped results might be not consistent with gene order
# THINK if multihitdb better
# FIX genes ids retrieval
# FIX to not get clusters from different genomes
# CHECK all the scores
def find_clusters(mapped_res, old_query_upd_scores, d_strand_flip_penalty, s_0, score_x_i):

    print('finding clusters')
    mapped_results = mapped_res.res_map_to_header
    index_list = mapped_res.ind_list
    target_db_lookup = mapped_res.target_db_lookup
    target_db_h = mapped_res.target_db_h
    logging.debug(f"s_0 in find_cluster is {s_0}")
    # Be careful, for whatever reason python adds 1 blank row to _h, only lookup is safe for the length
    L_space = len(target_db_lookup.iloc[:, 1])

    # Here i should get a dict of scores for all the target proteins
    # First i make dict of targets from the results and corresponding queries. It should work as every target should only happen to be once in the mappes results, best match
    matches_t_q_ids_dict = pd.Series(mapped_results.query_ID.values,index=mapped_results.ID).to_dict()
    # Okay, i need actually the dict for all the target proteins where I would have '' if they have no match
    all_t_q_ids_dict = {n:"''" for n in target_db_h["ID"].values}
    all_t_q_ids_dict.update(matches_t_q_ids_dict)
    #print(all_t_q_ids_dict)
    # Now, let's have another dict where it would be scores instead of just queries for the target
    #old_query_upd_scores = {'QIN93248.1': 1, 'QIN93249.1': 2, 'QIN93250.1': 3, 'QIN93251.1': 4}
    s_0 = -1
    matches_t_scores_dict = {t:old_query_upd_scores[matches_t_q_ids_dict[t]] for t in matches_t_q_ids_dict.keys()}
    logging.debug(f"old_query_upd_scores {old_query_upd_scores}")
    #print(matches_t_scores_dict)
    # Now, i need to make a bigger dictionary including s_0 scores for the target proteins with no match
    all_t_score_dict = {n:s_0 for n in target_db_h["ID"].values}
    all_t_score_dict.update(matches_t_scores_dict)
    # Let's make dict with i corresponding to certain target prots
    all_i_t_dict = pd.Series(target_db_h.ID.values,index=target_db_h.index).to_dict()
    #print(all_i_t_dict)

    
    # Here I should encode strand flips in the dict
    target_db_h_strand_list = target_db_h["strand"].values.tolist()
    strand_flips_dict = {k:1 if target_db_h_strand_list[k]!=target_db_h_strand_list[k-1] else 0 for k in range(0, L_space)}
    # just to make sure the first gene has 0 for strand flip
    strand_flips_dict[0] = 0
    #print('strand_flips_dict', strand_flips_dict)
    
    # I should also encode changes of genome in the dict
    # that makes genome id out of prot id, e.g. NC_111.1_1 -> NC_111.1
    genomes_col = target_db_h['ID'].str.rsplit('_', 1).str.get(0).tolist()
    genome_change_dict = {k:1 if genomes_col[k]!=genomes_col[k-1] else 0 for k in range(0, L_space)}

    cluster_matches_i_0 = list()
    cluster_matches_i_1 = list()
    cluster_matches_s_max = list()

    if d_strand_flip_penalty is None:
        d_strand_flip_penalty = 1
    else:
        d_strand_flip_penalty = -1*d_strand_flip_penalty
    # CHECK if score max cluster set up correct
    # Starting the Johannes algorithm here
    S_max = 0
    i_0 = 0
    i_1 = 0
    # this was not defined in latex
    S_min = 0
    S_i_minus_1 = 0
    S_i = 0
    i = 0
    # weird parameter not as in latex
    if_prev_clu_complete = 0

    while i < L_space:
        target_prot_id = all_i_t_dict[i]
        s_x_i = all_t_score_dict[target_prot_id]
        print('target_prot_id, s_x_i', target_prot_id, s_x_i)
        f_strand_flip = strand_flips_dict[i]
        print('S_i_minus_1, f_strand_flip, d_strand_flip_penalty', S_i_minus_1, f_strand_flip, d_strand_flip_penalty)
        genome_changed = genome_change_dict[i]
        if genome_changed == 1:
            S_max = 0
            S_min = 0
            S_i_minus_1 = 0
            S_i = 0

        if S_i_minus_1 - f_strand_flip*d_strand_flip_penalty + s_x_i > max(0, s_x_i):
            print('first')
            S_i = S_i_minus_1 - f_strand_flip*d_strand_flip_penalty + s_x_i
            print('S_i, S_max', S_i, S_max)
            if S_i > S_max:
                print('first second if')
                S_max = S_i
                i_1 = i
                print('i_1 = i', i_1, i)
            # not as in latex
            #else:
            #    S_max = s_x_i
            # end
        else:
            print('second')
            print('S_max, S_min', S_max, S_min)
            if S_max > S_min:
                print('second, append')
                cluster_matches_i_0.append(i_0)
                cluster_matches_i_1.append(i_1)
                cluster_matches_s_max.append(S_max)
                S_max = 0
                print(cluster_matches_i_0)
                print(cluster_matches_i_1)
            while s_x_i <= 0 and i < L_space - 1:
                #i = i + 1
                target_prot_id = all_i_t_dict[i]
                s_x_i = all_t_score_dict[target_prot_id]
                print('while, s_x_i', s_x_i, target_prot_id)
                genome_changed = genome_change_dict[i]
                if genome_changed == 1:
                    print('genome_changed')
                    S_max = 0
                    S_min = 0
                    S_i_minus_1 = 0
                    S_i = 0
                #    i = i - 1
                    break
                if s_x_i > 0:
                    i = i -1 
                    s_x_i = 0
                    break
                i = i + 1
            # not as in latex
            #if S_max != 0:
            #S_max = s_x_i
            # end of not as in latex
            S_i = s_x_i
            i_0 = i
            print('i_0', i_0)
        # not as in latex
        #if S_i <= 0:
        #    print('S_i <= 0')
        #    S_i = s_x_i
        #S_max = s_x_i
        S_i_minus_1 = S_i
        print('end S_i_minus_1', S_i_minus_1)
        i = i + 1
    if S_max > S_min:
        cluster_matches_i_0.append(i_0)
        cluster_matches_i_1.append(i_1)
        cluster_matches_s_max.append(S_max)

    # Now let's get the actual info for the clusters
    # MAKE FASTER
    print('cluster_matches_i_0', cluster_matches_i_0, 'cluster_matches_i_1', cluster_matches_i_1)
    targets1 = [all_i_t_dict[i] for i in cluster_matches_i_0]
    targets2 = [all_i_t_dict[i] for i in cluster_matches_i_1]
    print('targets1', targets1)
    print('targets2', targets2)
    #target_db_h_coord1_list = target_db_h["coord1"].values.tolist()
    #target_db_h_coord2_list = target_db_h["coord2"].values.tolist()
    cluster_starts = target_db_h["coord1"].iloc[cluster_matches_i_0].astype('int').values.tolist()
    cluster_ends = target_db_h["coord2"].iloc[cluster_matches_i_1].astype('int').values.tolist()
    #print(cluster_starts)
    #print('cluster_ends', cluster_ends)
    # Okay, now I have to make extended lists because i need all members of the cluster, not just left and right proteins
    cluster_matches_pd_ids = [[y for y in range(cluster_matches_i_0[w],cluster_matches_i_1[w]+1)] for w in range(0,len(cluster_matches_i_0))]
    #print(cluster_matches_pd_ids)
    #target_genes_ids = [[all_i_t_dict[y] for y in range(cluster_matches_i_0[w],cluster_matches_i_1[w]+1)] for w in range(0,len(cluster_matches_i_0))]
    # Now, let's better make one for loop to make all the lists with info
    target_genes_ids = []
    prots_strands = []
    coords1 = []
    coords2 = []
    query_genes_ids = []
    cluster_matches = []
    for w in range(0, len(cluster_matches_pd_ids)):
        current_ids_set = cluster_matches_pd_ids[w]
        target_genes_ids.append([all_i_t_dict[y] for y in cluster_matches_pd_ids[w]])
        prots_strands.append(target_db_h["strand"].iloc[cluster_matches_pd_ids[w]].astype('int').values.tolist())
        coords1.append(target_db_h["coord1"].iloc[current_ids_set].astype('int').values.tolist())
        coords2.append(target_db_h["coord2"].iloc[current_ids_set].astype('int').values.tolist())
        #print(target_genes_ids[w])
        #print([all_t_q_ids_dict[s] for s in target_genes_ids[w]])
        query_genes_ids.append([all_t_q_ids_dict[s] for s in target_genes_ids[w]])
        #cluster_matches.append((i_0_cluster_start, i_1_cluster_end, score_max_cluster, query_genes_ids, target_genes_ids, prots_strands, coords1, coords2))
    cluster_matches_df = pd.DataFrame({0:cluster_starts, 1:cluster_ends, 2:cluster_matches_s_max, 3:query_genes_ids, 4:target_genes_ids, 5:prots_strands, 6:coords1, 7:coords2})
    print(cluster_matches_df)
    logging.debug(f"cluster_matches_df: \n {cluster_matches_df}")

    print(x)

    
    print(cluster_matches)
    print(len(cluster_matches))
    return cluster_matches


# it re-defines the scores which I have got in the first iteration 
# ADD significant clusters determination?
def update_scores_for_cluster_matches(significant_clusters_eval_filter_df, mapped_res, bias):

    # CHANGE for really significant clusters, not as it is?
    #significant_clusters = cluster_matches 
    # ADD query id to mapped results
    mapped_results = mapped_res.res_map_to_header

    target_db_lookup = mapped_res.target_db_lookup
    target_db_h = mapped_res.target_db_h

    # CHECK if these are right columns
    # CHECK if the query and target assignment is correct
    K = len(significant_clusters_eval_filter_df)
    L = len(target_db_lookup.index)

    logging.debug(f"L: {L}")
    
    # Should cluster prots be done better?
    cluster_prots = pd.DataFrame()
    cluster_prots['query_id'] = significant_clusters_eval_filter_df['query_prots'].explode()
    cluster_prots['target_id'] = significant_clusters_eval_filter_df['target_prots'].explode()

    l_orig = len(cluster_prots)
    logging.debug(f"l1: {l_orig}")
    # adding new pseudocounts here
    m_pseudocount = 2
    # each cluster should have number of prots + pseudocount
    l = l_orig + m_pseudocount*K
    logging.debug(f"l2: {l}")


    #logging.debug(f"K, L, l: {K, L, l}")
    
    #logging.debug(f"sign_clusters_df: \n {sign_clusters_df}")
    #logging.debug(f"cluster_prots: \n {cluster_prots}")
    #logging.debug(f"mapped results: \n {mapped_results}")
    # CHECK if it is correct to iterate through cluster matches?
    # So I am calculating enrichment score for cluster matches and for cluster prots with no match
    # CHECK if leaving only unique entries is correct
    # pseudocounts are added with parameter aplha_pseudocount
    aplha_pseudocount = pow(10,np.log10(1/L)-1)
    x_number_of_queries = len(mapped_results['query_ID'].unique())
    #logging.debug(f"x_number_of_queries: {x_number_of_queries}")
    # Do I need unique?
    for query_id in cluster_prots['query_id'].unique():
        logging.debug(f"query_id: {query_id}")
        # to avoid empty queries corresponding to prots with no match
        if query_id == '':
            significant_clusters_eval_filter_df.loc[significant_clusters_eval_filter_df['queries_string'].str.contains(str(query_id)),
        'list_new_scoreS_enrich'] = significant_clusters_eval_filter_df.loc[significant_clusters_eval_filter_df['queries_string'].str.contains(str(query_id)),
        'list_new_scoreS_enrich'].astype(str) + ',' + ''
            continue
        M_x = mapped_results['query_ID'][mapped_results['query_ID'] == query_id].shape[0]
        # There was a need to use this file (mapped_results_mish) here as the mapped results
        # have only best hits, not like initial m8. Also query there is in col 9
        #cmd = "cat mapped_results_mish | awk -F'\t' '{print $9}' | grep " + query_id +" | wc -l"
        #pr = subprocess.Popen(cmd, shell = True, stdout = subprocess.PIPE)
        #out, er = pr.communicate()
        #M_x = int(out.decode('utf-8'))
        m_x = cluster_prots[cluster_prots['query_id'] == query_id]['query_id'].count()
        # in this case M_x and m_x are equal as there were no single hits. 
        # CHECK with other data where would be hits not only in clusters
        logging.debug(f"M_x, m_x: {M_x, m_x}")
        # using pseudocounts for m_x/l proportion to avoid zeros in log
        cluster_prot_proportion = np.divide(m_x, l)
        cluster_prot_proportion = np.divide((m_x+aplha_pseudocount), (l + x_number_of_queries*aplha_pseudocount))
        score_x = np.log(np.divide(cluster_prot_proportion, np.divide(M_x, L))) - bias
        logging.debug(f"updated score for q: {query_id} is {score_x}")


        old_query_upd_scores[query_id] = score_x

        for target_id in cluster_prots['target_id'][cluster_prots['query_id'] == query_id]:
            if iter_counter == 1:
                significant_clusters_eval_filter_df.loc[significant_clusters_eval_filter_df['targets_string'].str.contains(str(target_id)),
                    'initial_q_or_match'] = True
            else:
                if query_id in q_and_matches:
                    significant_clusters_eval_filter_df.loc[significant_clusters_eval_filter_df['targets_string'].str.contains(str(target_id)),
                    'initial_q_or_match'] = True

        # MAKE faster?
        # adding score of the query prot to get summarized score for the cluster
        # CHECK if correct
        significant_clusters_eval_filter_df.loc[significant_clusters_eval_filter_df['queries_string'].str.contains(str(query_id)),
        'new_score_enrich'] = significant_clusters_eval_filter_df.loc[significant_clusters_eval_filter_df['queries_string'].str.contains(str(query_id)),
        'new_score_enrich'] + score_x
        significant_clusters_eval_filter_df.loc[significant_clusters_eval_filter_df['queries_string'].str.contains(str(query_id)),
        'list_new_scoreS_enrich'] = significant_clusters_eval_filter_df.loc[significant_clusters_eval_filter_df['queries_string'].str.contains(str(query_id)),
        'list_new_scoreS_enrich'].astype(str) + ',' + str(score_x)

        #print(sign_clusters_df['list_new_scoreS_enrich'])
        # FIX figure out how to set b and threshold for prot to be enriched in cluster
        # ADD file to store scores corresponding to profiles
    # calculating the scores for cluster prots with no match
    # THINK if it is ok to iterate again and also the order of scores adding to the
    # list of scores is not competely correct
    matches_ids_list = mapped_results['ID'].tolist()

    # Counting all target cluster prots having query matches. The first line had problems as it counted all the proteins in the clusters, also the one that had no matches to the query
    #m_x_sum = len(cluster_prots['query_id'])
    non_match_counter = len(cluster_prots[cluster_prots['query_id'] == ''])  
    m_x_sum = l_orig - non_match_counter
    logging.debug(f"m_x_sum is: {m_x_sum} , non_match_counter is {non_match_counter}")
    # Counting all target proteins with matches to some queries
    M_x_sum = len(mapped_results['ID'])
    #s_0 = np.log(np.divide(np.divide((l-m_x_sum),l), np.divide((L-M_x_sum),L))) - bias
    #print(x)
    # added pseudocount for the rare case when there are no gaps in cluster
    # Other pseudocounts are added to have adequate value when all prots in cluster are matches
    s_0 = np.log(np.divide(np.divide((l-m_x_sum+aplha_pseudocount),l + x_number_of_queries*aplha_pseudocount), np.divide((L-M_x_sum),L))+aplha_pseudocount) - bias
    logging.debug(f's_0 {s_0}, m_x_sum {m_x_sum}, M_x_sum {M_x_sum}, L {L}, l {l} ')
    print(f's_0 {s_0}, m_x_sum {m_x_sum}, M_x_sum {M_x_sum}, L {L}, l {l} ')
    for target_id in cluster_prots['target_id']:
        logging.debug(f"the cycle is running")
        if target_id not in matches_ids_list:
            logging.debug(f"not in list")
            logging.debug(f"no match: {target_id}")
            #s_0 = np.divide(np.divide((l-m_x_sum),l), np.divide((L-M_x_sum),L)) - bias
            logging.debug(f"updates s_0 for: {target_id} is= {s_0}")
            #logging.debug(f"sign_clusters_df[pd.DataFrame(sign_clusters_df['target_prots'].tolist()).isin(target_id.split()).any(1).values]")
            #tmp_df = sign_clusters_df[pd.DataFrame(sign_clusters_df['target_prots'].tolist()).isin(target_id.split()).any(1).values]['new_score_enrich'] + s_0
            array_of_bool = pd.DataFrame(significant_clusters_eval_filter_df['target_prots'].tolist()).isin(target_id.split()).any(axis=1).values
            index_of_row = int(np.where(array_of_bool == True)[0])
            significant_clusters_eval_filter_df.loc[significant_clusters_eval_filter_df.index[index_of_row],'new_score_enrich'] = significant_clusters_eval_filter_df.loc[significant_clusters_eval_filter_df.index[index_of_row],'new_score_enrich'] + s_0
            # FIX to make the cell string 
            #sign_clusters_df.loc[sign_clusters_df.index[index_of_row],'list_new_scoreS_enrich'] = sign_clusters_df.loc[sign_clusters_df.index[index_of_row],'list_new_scoreS_enrich'] + ',' + str(s_0)
            #logging.debug(f"sign_clusters_df: \n {sign_clusters_df}")
            #logging.debug(f"sign_clusters_df['list_new_scoreS_enrich'].tolist(): {sign_clusters_df['list_new_scoreS_enrich'].tolist()}")
        # I have NO IDEA what these lines here are supposed to be for
        #else:
        #    logging.debug(f"is in list")
        #    s_0 = -1 - bias
    #if len(cluster_prots['target_id']) == 0:
    #    logging.debug(f"0 target in clusters")
    #    s_0 = -1 - bias
    logging.debug(f"s_0: {s_0}")
    print('significant_clusters_eval_filter_df2 \n', significant_clusters_eval_filter_df)
    logging.debug(f"sign_clusters_df: \n {significant_clusters_eval_filter_df}")
    significant_clusters_eval_filter_df.to_csv(files.res + '_' + str(iter_counter) + '_sign_clusters_df_filter_upd', sep = '\t', index = False)
    return(significant_clusters_eval_filter_df, s_0, old_query_upd_scores, L, l)


# CHANGE to not duplicate so much the func above
def calculate_karlin_stat(cluster_matches, mapped_res, s_0, bias):
    print('using c code for Karlin-Altschul statistics')

    curr_file_path = os.path.realpath(__file__)
    logging.debug(f"curr_file_path: {curr_file_path}")
    karlin_file_path = '/'.join(curr_file_path.split('/')[:-1]) + '/karlin_c.c'
    subprocess.call(['cc', '-fPIC', '-shared', '-o', 'karlin_c.so', karlin_file_path])
    # !! CHANGE the path later
    # THINK to do it only once
    so_file = "./karlin_c.so"
    my_functions = CDLL(so_file)

    # ASK Johannes if that's ok that I calculated enrich scores for all search results,
    # not only the clusters

    mapped_results = mapped_res.res_map_to_header

    logging.debug(f"calculating prob for each match")
    logging.debug(f"mapped_results: \n {mapped_results}")

    significant_clusters = cluster_matches 

    # CHECK if these are right columns
    # CHECK if the query and target assignment is correct
    K = len(significant_clusters)
    
    target_db_lookup = mapped_res.target_db_lookup
    L = len(target_db_lookup.index)
    
    sign_clusters_df = pd.DataFrame(significant_clusters)
    sign_clusters_df.columns = ["coord1", "coord2", "score",
    "query_prots", "target_prots", "strand", "coordS1", "coordS2"]

    cluster_prots = pd.DataFrame()
    cluster_prots['query_id'] = sign_clusters_df['query_prots'].explode()
    cluster_prots['target_id'] = sign_clusters_df['target_prots'].explode()

    l = len(cluster_prots)
    logging.debug(f"len(cluster_prots): {l}")
    aplha_pseudocount = pow(10,np.log10(1/L)-1)

    # Now let's start to do the karlin-altschul statistics properly, not with the updated score but with the target prots scores that are derived from the corresponding query scores (not yet updated)
    all_target_prots = pd.DataFrame()
    all_target_prots['ID'] = target_db_lookup['ID']
    matched_target_prots = mapped_results[['ID', 'query_ID']]
    matched_target_prots['score_x'] = matched_target_prots['query_ID'].map(old_query_upd_scores)
    
    all_target_prots['query_ID'] = 0
    # That is to artificially round s_0 not to 0 if s_0 is negative (otherwise it violates karlin-altschul requiremens of the expected negative score)
    if s_0 > -0.5 and s_0 < 0:
        s_0 = -0.51
    all_target_prots['score_x'] = s_0
    all_target_prots.update(matched_target_prots)
    enrich_scores = all_target_prots['score_x'].to_numpy()
    print(s_0)
    print(enrich_scores)
    
    logging.debug(f"scores for all res: \n {enrich_scores}")
    logging.debug(f"len: {len(enrich_scores)}")

    #if 0 not in enrich_scores:
    #    no_0 = True
    #    enrich_scores = np.append(enrich_scores, 0)
    #else:
    #    no_0 = False
    #print('no_0', no_0)

    # REMOVE later, just for current test

    logging.debug(f"enrich_scores: \n {enrich_scores}")
    # ASK Johannes if I done the scores correctly
     # # expand unique scores in order to have larger range of integer, wider range to calc e-value
    # # multipliying param for log-scores
    #mult_param = 6
    # THINK if I should give an option to set the multiplying parameter to the user
    mult_param = 1
    scores_bins = np.histogram_bin_edges(enrich_scores, bins='auto')
    logging.debug(f"BINS: {scores_bins}")
    logging.debug(f"bins diff: {scores_bins[1]-scores_bins[0]}")
    # # THINK!

    # ASK Johannes, REMOVE later, that's just to fix problem with log when m_x = 0
    # ASK Johannes if it is ok to use binning, no I think, it is NOT ok
    #unique_scores = mult_param*np.unique(enrich_scores)

    sorted_scores = np.sort(enrich_scores)
    logging.debug(f"multiplied, sorted: \n {sorted_scores}")

    sorted_scores = np.round(sorted_scores, decimals=0)
    unique_scores, score_counts = np.unique(sorted_scores, return_counts=True)

    unique_scores = unique_scores.astype(int)
    logging.debug(f"unique_scores: \n {unique_scores}")
    logging.debug(f"uniq enrich scores: \n {np.unique(enrich_scores)}")
    logging.debug(f"{np.unique(mult_param*enrich_scores)}")
    logging.debug(f"{np.round(unique_scores, decimals=0)}")
    logging.debug(f"len of orig scores {len(np.unique(enrich_scores))}")
    logging.debug(f"len of scores {len(unique_scores)}")

    logging.debug(f"score_counts {score_counts}")
    logging.debug(f"len {len(enrich_scores)}")

    #score_prob = score_counts / len(enrich_scores)
    # In fact, here should be all the proteins from all the genomes (incl those without matches)
    #score_counts = score_counts + (target_ids_counts - 1)
    score_prob = score_counts / L
    logging.debug(f"score_prob {score_prob}")

    # figure out why round works weird so -3.5 > 4, -2.5 > 2
    scores_table = np.column_stack((np.round(unique_scores, decimals=0), score_prob))
    # sort so that values would go in correct order for inserting the missing ones
    scores_table = scores_table[scores_table[:, 0].argsort()]
    logging.debug(f"scores_table \n {scores_table}")
    # MAKE faster!
    # the loop to insert missed integers (should be one by one) to scores and 0 to correspond. probs
    prev_int_val = int(scores_table[0][0])
    curr_row_ind = 0
    # CHECK if the table is sorted
    for i in scores_table:
        logging.debug(f"start inserting")
        logging.debug(f"int_val {prev_int_val}")
        logging.debug(f"i {i}")
        # CHECK what negative value in times_insert would do (seems like does nothing correctly)
        times_insert = int(i[0]) - prev_int_val - 1
        logging.debug(f"times_insert {times_insert}")
        
        for j in range(0, times_insert):
            insert_index = np.where(np.all(scores_table==i,axis=1))[0][0]
            logging.debug(f"insert_index {insert_index}")
            int_to_insert = prev_int_val + j + 1
            scores_table = np.insert(scores_table, insert_index, np.array((int_to_insert, 0)), 0)  
            j = j + 1
        prev_int_val = int(i[0])
        curr_row_ind = curr_row_ind + 1
        logging.debug(f"scores_table {scores_table}")

    # REMOVE!! tmp! for little dataset testing to have 0
    #scores_table = np.insert(scores_table, 1, np.array((0, 0.33)), 0)  
    #scores_table = np.insert(scores_table, 2, np.array((1, 0.17)), 0)  
    #print(scores_table)
    
    unique_scores = scores_table[:,0]
    score_prob = scores_table[:,1]

    # ASK Johannes what to do if all the scores are negative?
    print(unique_scores)
    index_of_0 = np.where(unique_scores == 0)[0][0]
    
    logging.debug(f"final uniq and prob \n {unique_scores} \n {score_prob}")
    arr = score_prob
    
    # ASK Johannes why I need log2 and where to make these scores conversions
    # THINK if it is slow
    # FIX how it works, all the scores should be in this format, no just ignored
    # THINK if I need it given I have it above, maybe remove it

    # if unique_scores[0] != 0:
    #     if unique_scores[0] >0:
    #         min_score = round(array_width_koef*np.log2(unique_scores[0]))
    #     else:
    #         min_score = round(array_width_koef*np.log2(abs(unique_scores[0])))*-1
    # else:
    #     min_score = 0
    # max_score = round(array_width_koef*np.log2(unique_scores[len(unique_scores)-1]))

    # REMOVE or FIX, all scores should be in the same format (log and so on)
    min_score = int(np.min(unique_scores))
    max_score = int(np.max(unique_scores))
    logging.debug(f"min and max {min_score, max_score}")
    
    logging.debug(f"index_of_0 {index_of_0}")
    logging.debug(f"prob of 0 {arr[index_of_0:]}")
    logging.debug(f"arr {arr}")
    arr = np.array(arr)
    logging.debug(f" {type(arr)}")
    logging.debug(f"element type {type(arr[0])}")
    
    pointer_to_middle = arr[index_of_0:].ctypes.data_as(ct.POINTER(ct.c_double))
    pointer_to_middle = arr[index_of_0:].ctypes.data_as(ct.POINTER(ct.c_double))
    logging.debug(f"{pointer_to_middle.contents}")

    
    my_functions.BlastKarlinLambdaNR.argtypes = [POINTER(c_double),c_int32, c_int32]
    my_functions.BlastKarlinLambdaNR.restype = c_double
    BlastKarlin_lambda = my_functions.BlastKarlinLambdaNR(pointer_to_middle, min_score, max_score)
    logging.debug(f"lambda {BlastKarlin_lambda}")
    BlastKarlin_lambda = ct.c_double(BlastKarlin_lambda)
    
    my_functions.BlastKarlinLtoH.restype = c_double
    BlastKarlin_H = my_functions.BlastKarlinLtoH(pointer_to_middle, min_score, max_score,
     BlastKarlin_lambda)
    logging.debug(f"H {BlastKarlin_H}")
    BlastKarlin_H = ct.c_double(BlastKarlin_H)
    

    my_functions.BlastKarlinLHtoK.restype = c_double
    BlastKarlin_K = my_functions.BlastKarlinLHtoK(pointer_to_middle, min_score, max_score,
     BlastKarlin_lambda, BlastKarlin_H)
    logging.debug(f"K {BlastKarlin_K}")
    return(BlastKarlin_lambda, BlastKarlin_K)



# ASK where strand flip penalty goes? next iter?
def set_strand_flip_penalty(cluster_matches, mapped_res):
    target_db_h = mapped_res.target_db_h
    mapped_results = mapped_res.res_map_to_header
    #print(target_db_h)
    # this just delete the consecutive duplicates, therefore I can cound 
    # how many switches between strands are in the file (1 1 -1 -1 1 -> 1 -1 1)
    F = len([key for key, _group in itertools.groupby(target_db_h["strand"])]) - 1
    logging.debug(f"F = {F}")
    # CHANGE to not duplicate update_scores
    significant_clusters = cluster_matches
    sign_clusters_df = pd.DataFrame(significant_clusters)
    sign_clusters_df.columns = ["coord1", "coord2", "score",
    "query_prots", "target_prots", "strand", "coordS1", "coordS2"]
    logging.debug(f"sign_clusters_df \n {sign_clusters_df}")
    f = 0
    # this should count number of flips in each cluster
    # think how to speed up this step
    for cluster_strands in sign_clusters_df["strand"]:
        n_flips_cluster = len([key for key, _group in itertools.groupby(cluster_strands)]) - 1
        f = f + n_flips_cluster
    logging.debug(f"f = {f}")
    # CHANGE to not duplicate scores updating
    K = len(significant_clusters)
    cluster_prots = pd.DataFrame()
    cluster_prots['query_id'] = sign_clusters_df['query_prots'].explode()
    l = len(cluster_prots)
    print("K, l = ", K, l)
    target_db_lookup = mapped_res.target_db_lookup
    L = len(target_db_lookup.index)

    # adding new pseudocounts here to allow more flips
    strand_flip_pseudocount = 1
    # each cluster should have number of prots + pseudocount
    f_corrected = f + strand_flip_pseudocount
    # ASK Johannes if pseudocount is correct
    aplha_pseudocount = pow(10,np.log10(1/L)-1)
    # CHECK if it is right that in case of all res probs I should iterate through target
    # NOT through query as in example above
    n_target = len(mapped_results['ID'].unique())

    #if f == 0:
    #    f = F/100000
    cluster_flips_proportion = np.divide(f_corrected, (l - K))
    cluster_flips_proportion = np.divide(f_corrected + aplha_pseudocount, (l + n_target*aplha_pseudocount - K))
    cluster_flips_proportion = np.divide(f_corrected + aplha_pseudocount, l)
    d = np.log(np.divide(cluster_flips_proportion, np.divide(F, (L - 1))))
    logging.debug(f'd = {d} f, f_corrected, F, l, K, L {f, f_corrected, F, l, K, L}')
    # is it a good place to return?
    # CHECK if these are really significant clusters
    return(d)


# ASK Johannes if in e-value calculations the L should actually be the L, not like L*L
def calculate_e_value(stat_lambda, stat_K, cluster_matches, mapped_res):
    # FIX to be variable taken from number of prots in target
    # CHANGE it to be linked and the same with the used in calculate_karlin_stat
    # THINK if I should give an option to set the multiplying param to an user
    mult_param = 1
    print('calculating e-value')
    target_db_lookup = mapped_res.target_db_lookup
    L = len(target_db_lookup.index)


    sign_clusters_df = pd.DataFrame(cluster_matches)
    sign_clusters_df.columns = ["coord1", "coord2", "score",
    "query_prots", "target_prots", "strand", "coordS1", "coordS2"]
    sign_clusters_df["initial_q_or_match"] = False

    logging.debug(f"sign_clusters_df: \n {sign_clusters_df}")
    sign_clusters_df['new_score_enrich'] = 0
    sign_clusters_df['list_new_scoreS_enrich'] = ''

    # MAKE faster?
    sign_clusters_df['queries_string'] = [','.join(map(str, l)) for l in sign_clusters_df['query_prots']]
    sign_clusters_df['targets_string'] = [','.join(map(str, l)) for l in sign_clusters_df['target_prots']]
    logging.debug(f"sign_clusters_df \n {sign_clusters_df}")
    print(sign_clusters_df)
    
    stat_lambda = stat_lambda.value
    logging.debug(f" {stat_lambda, stat_K}")
    sign_clusters_df["e-value"] = 0
    # PUT in correct place
    pd.options.display.float_format = '{:,.4f}'.format
    for r in range(0, len(sign_clusters_df.index)):
        logging.debug(f"r {r}")
        enrich_score = sign_clusters_df["score"][r]
        # converting enrich score as for karlin stats calculation
        conv_enrich_score = int(np.round(enrich_score*mult_param, decimals = 0))
        logging.debug(f"conv_enrich_score {conv_enrich_score}")
        # is it okay if that's a sum?
        evalue = stat_K*L*np.exp(stat_lambda*(-1)*conv_enrich_score)
        sign_clusters_df.iat[r, sign_clusters_df.columns.get_loc("e-value")] = float(evalue)
        logging.debug(f"eval = {evalue}")
    # here I will try to calculate the q-value suggested by Johannes and replace e-value with them. It is a temporary option, just to try out
    # ----disabled for now as I have to figure out the problem with e-val
    #sign_clusters_df_sorted_by_eval = sign_clusters_df.sort_values('e-value')
    #sign_clusters_df_sorted_by_eval['divided'] = sign_clusters_df_sorted_by_eval['e-value'].divide(sign_clusters_df.index + 1, axis="rows")
    #sign_clusters_df_sorted_by_eval['min_fdr'] = np.minimum.accumulate(sign_clusters_df_sorted_by_eval.divided[::-1])
    #print('sign_clusters_df_sorted_by_eval', sign_clusters_df_sorted_by_eval)
    #sign_clusters_df_sorted_by_eval['e-value'] = sign_clusters_df_sorted_by_eval['min_fdr']
    #sign_clusters_df = sign_clusters_df_sorted_by_eval.drop(columns=['divided', 'min_fdr'])
    #print(sign_clusters_df)
    # -----------
    # CHANGE the name? 
    # they are not filtered by e-val yet
    cluster_path = files.res + '_' + str(iter_counter) + '_iter_sign_clusters_enrich_stat'
    sign_clusters_df.to_csv(cluster_path, sep = '\t', index = False)
    # Check if e-val = 0.01 is good filtering
    eval_filter = float(args.eval)
    use_eval_filter = args.evalfilteruse
    if use_eval_filter == '1':
        significant_clusters_eval_filter_df = sign_clusters_df.loc[(sign_clusters_df['e-value'] <eval_filter) & (sign_clusters_df['e-value'] >= 0)]
        cluster_path2 = files.res + '_' + str(iter_counter) + '_iter_sign_clusters_enrich_stat_filtered'
        significant_clusters_eval_filter_df.to_csv(cluster_path2, sep = '\t', index = False)
    else:
        print('e-value filter disabled')
        significant_clusters_eval_filter_df = sign_clusters_df.copy()
    logging.debug(f"significant_cluster_df_enriched \n {sign_clusters_df}")
    logging.debug(f"significant_clusters_eval_filter_df \n {significant_clusters_eval_filter_df}")
    return sign_clusters_df, significant_clusters_eval_filter_df


def extract_proteins_cluster_neighborhood(sign_clusters_df, mapped_res):
    # that extracts not only the neighbourhood but also what is inside the cluster
    # here i refer to neighbourhood often also including the non-matches from inside the clusters
    print('updating query profile started (extracting prots)')
    logging.debug(f"significant (?) clusters table \n {sign_clusters_df}")
    target_fasta = args.targetfa
    logging.debug(f"{sign_clusters_df['target_prots']}")
    neighbourhood_path = files.res + '_' + str(iter_counter) + '_iter_target_clusters_neighbourhood'

    # extracting 3 neghbours from each side
    # MAKE a parameter
    neighbors_number_1side = 3
    if if_singleton == '1':
        neighbors_number_1side = 0

    # WHY do I need this file?
    #target_protID_cluster_file_idx = open("target_protID_cluster_file_idx", "w")
    target_db_lookup = mapped_res.target_db_lookup
    #print(target_db_lookup)
    print('creating subdb')

    # arrays with the left and right proteins of each cluster
    arr_prot_id_left = sign_clusters_df["target_prots"].str[0].to_numpy()
    arr_prot_id_right = sign_clusters_df["target_prots"].str[-1].to_numpy()
    arr_entire_cluster_id = sign_clusters_df["target_prots"].explode('target_prots').to_numpy()
    # array with the protein real ids taken from target db lookup. Indices in array can be
    # used to access the mmseqs ids (ind_id)
    arr_target_db_lookup_real_ids = target_db_lookup['ID'].to_numpy()
    arr_target_db_lookup_mmseqs_ind = target_db_lookup['ind_id'].to_numpy()

    # i need to separate matches and the rest of the cluster + neighbourhood
    matches_ids = mapped_res.res_map_to_header['ID'].to_numpy()
    arr_matches_in_clu = np.intersect1d(arr_entire_cluster_id,matches_ids)

    # make array with the neighbours indices
    l_all_indices_clu_neigh = []
    clu_indices_for_frac_occ_min = []
    init_clu_ind = 0
    target_lookup_len = len(arr_target_db_lookup_real_ids)
    print('target_lookup_len', target_lookup_len)
    for i in range(arr_prot_id_left.size):
        left_ind = int(np.where(arr_target_db_lookup_real_ids == arr_prot_id_left[i])[0])
        right_ind = int(np.where(arr_target_db_lookup_real_ids == arr_prot_id_right[i])[0])
        left_border = left_ind-neighbors_number_1side
        right_border = right_ind+neighbors_number_1side+1
        # that is for the edge cases, protein id must be prodigal_derived to get the genome id
        genome_id = '_'.join(arr_prot_id_left[i].split('_')[:-1])
        last_ind_for_this_genome = target_db_lookup[target_db_lookup['ID'].str.contains(genome_id)].index[-1]
        first_ind_for_this_genome = target_db_lookup[target_db_lookup['ID'].str.contains(genome_id)].index[0]
        #print('arr_prot_id_left[i]', arr_prot_id_left[i], 'arr_prot_id_right[i]', arr_prot_id_right[i], 'last_ind_for_this_genome', last_ind_for_this_genome, 'first_ind_for_this_genome', first_ind_for_this_genome)
        #print('left_border right_border', left_border, right_border)
        if left_border < first_ind_for_this_genome:
            l_all_indices_clu_neigh.extend(list(range(first_ind_for_this_genome,right_border+1)))
            l_all_indices_clu_neigh.extend(list(range((last_ind_for_this_genome-(first_ind_for_this_genome - left_border)+1),last_ind_for_this_genome+1)))
        elif right_border > last_ind_for_this_genome:
            l_all_indices_clu_neigh.extend(list(range(left_border,last_ind_for_this_genome+1)))
            l_all_indices_clu_neigh.extend(list(range(first_ind_for_this_genome,first_ind_for_this_genome+(right_border-last_ind_for_this_genome))))
            #print('arr_prot_id_right', arr_prot_id_right[i], 'genome_id', genome_id, 'first_ind_for_this_genome', first_ind_for_this_genome)
            #print('l_all_indices_clu_neigh', l_all_indices_clu_neigh)
        # here is the normal expected in most cases situation, not the edge case
        else:
            l_all_indices_clu_neigh.extend(list(range(left_border,right_border+1)))
        clu_indices_for_frac_occ_min.extend([init_clu_ind]*len(range(left_border,right_border+1)))
        init_clu_ind = init_clu_ind + 1
        print('l_all_indices_clu_neigh', l_all_indices_clu_neigh)
        print('clu_indices_for_frac_occ_min', clu_indices_for_frac_occ_min)
    logging.debug(f"l_all_indices_clu_neigh {l_all_indices_clu_neigh}")
    # let's get the same indices but now for the matches in cluster only
    indices_matches_in_clu = np.unique(np.where(np.isin(arr_target_db_lookup_real_ids, arr_matches_in_clu)))
    # add indices of proteins between left and right prots of the cluster and keep uniq
    all_indices_clu_neigh = np.unique(np.array(l_all_indices_clu_neigh))
    # also for another function i would get the protein ids for these
    arr_clu_neigh_prots = arr_target_db_lookup_real_ids[all_indices_clu_neigh]
    # here I need a non-unique subset of the real proteins to correspond to clu_indices_for_frac_occ_min array, as I need for the clu indices a dataframe with the corresponding real ids to merge later with the seq clu results in add_new_profiles to see how many spatial clusters were for each sequences cluster
    arr_clu_neigh_prots_non_uniq = arr_target_db_lookup_real_ids[l_all_indices_clu_neigh]
    #print('l_all_indices_clu_neigh', l_all_indices_clu_neigh)
    clu_indices_for_frac_occ_min_df = pd.DataFrame({'real_prot_id': arr_clu_neigh_prots_non_uniq, 'internal_clu_ind': clu_indices_for_frac_occ_min})
    #print('clu_indices_for_frac_occ_min_df', clu_indices_for_frac_occ_min_df)
    # obtain mmseqs ids corresponding to indices
    arr_mmseqs_ind_clu_neigh = arr_target_db_lookup_mmseqs_ind[all_indices_clu_neigh]
    # these are indices for the matches in the cluster only
    arr_mmseqs_ind_matches_in_clu = arr_target_db_lookup_mmseqs_ind[indices_matches_in_clu]

    # now I have to finally remove matches indices from the rest
    arr_mmseqs_ind_clu_neigh_only = np.setdiff1d(arr_mmseqs_ind_clu_neigh,arr_mmseqs_ind_matches_in_clu)

    # I still want to have the file with all clusters + neighbourhood proteins as I need them later
    # However, I do not proceed here the same way with the matches/non-matches only files, i do it in keep_enriched as I need to filter them before in initialize_new_prot_score2. I will just pass the np arrays
    np.savetxt('target_protID_cluster_file_idx', arr_mmseqs_ind_clu_neigh, fmt='%i')
 
    #target_protID_cluster_file_idx.close()
    with open('target_protID_cluster_file_idx_sorted','w') as out0:
        subprocess.call(['sort', '-u', '-n', 'target_protID_cluster_file_idx'],stdout=out0)

    # BE careful! this db accessory files are not in order with db seqs
    subprocess.call(['mmseqs', 'createsubdb', 'target_protID_cluster_file_idx_sorted', 
     args.targetdb, str(neighbourhood_path)+'_db'])
    subprocess.call(['mmseqs', 'createsubdb', 'target_protID_cluster_file_idx_sorted', 
     str(args.targetdb)+'_h', str(neighbourhood_path)+'_db_h'])
    subprocess.call(['mmseqs', 'convert2fasta', str(neighbourhood_path)+'_db', 
     neighbourhood_path])

    # Here there are files with the mmseqs ids for only matches (to update profiles later) and only neighbourhood (to make new profiles)
    # Even though i call it neigh_only, these are not only neighbours but also the non-matches inside of the clusters.
    np.savetxt('target_protID_cluster_file_idx_neigh_only', arr_mmseqs_ind_clu_neigh_only, fmt='%i')
    np.savetxt('target_protID_cluster_file_idx_matches_only', arr_mmseqs_ind_matches_in_clu, fmt='%i')

    with open('target_protID_cluster_file_idx_neigh_only_sorted','w') as out:
        subprocess.call(['sort', '-u', '-n', 'target_protID_cluster_file_idx_neigh_only'],stdout=out)
    with open('target_protID_cluster_file_idx_matches_only_sorted','w') as out1:
        subprocess.call(['sort', '-u', '-n', 'target_protID_cluster_file_idx_matches_only'],stdout=out1)

    # I do not convert to fasta for these (clu matches/non-matches) as I dont need it later
    subprocess.call(['mmseqs', 'createsubdb', 'target_protID_cluster_file_idx_neigh_only_sorted', 
     args.targetdb, str(files.res) + '_' + str(iter_counter) +'_neigh_only_db'])
    subprocess.call(['mmseqs', 'createsubdb', 'target_protID_cluster_file_idx_neigh_only_sorted', 
     str(args.targetdb)+'_h', str(files.res) + '_' + str(iter_counter) +'_neigh_only_db_h'])
    
    subprocess.call(['mmseqs', 'createsubdb', 'target_protID_cluster_file_idx_matches_only_sorted', 
     args.targetdb, str(files.res) + '_' + str(iter_counter)+'_matches_only_db'])
    subprocess.call(['mmseqs', 'createsubdb', 'target_protID_cluster_file_idx_matches_only_sorted', 
     str(args.targetdb)+'_h', str(files.res) + '_' + str(iter_counter)+'_matches_only_db_h'])

    #target_clusters_neighbourhood.close()
    return arr_mmseqs_ind_matches_in_clu, arr_mmseqs_ind_clu_neigh_only, arr_clu_neigh_prots, arr_matches_in_clu, clu_indices_for_frac_occ_min_df


def reassign_non_enriched(old_query_upd_scores, bias, s_0):
    # the function is used to filter the query profiles (representatives of the clusters) by their enrichment
    # the proteins (profile representatives) that are under the threshold of enrichment/bias should be assigned s_0 as the score (as for the gaps)
    print('old_query_upd_scores', old_query_upd_scores)
    old_query_upd_scores = {k: s_0 if v < bias else v for (k, v) in old_query_upd_scores.items() }
    print('bias ', bias, ' s_0 ', s_0)
    print('old_query_upd_scores', old_query_upd_scores)



def update_profiles():
    # that is used if more than 1 iteration was run
    # here i will update the old msas with their matches and construct new msas with the proteins from the neighbourhood and make them to the one profile

    # Let's start with updating the old profiles
    #subprocess.call(['mmseqs', 'search', 
    #files.query_db + '_clu' + '_rep' + '_profile',
    #files.query_db + '_clu_msa_db_profile',
     #str(files.res) + '_' + str(iter_counter)+'_matches_only_db',
     #str(files.res) + '_' + str(iter_counter) + '_matches_only_db' + '_upd_res',
     #'tmp', '-a', '--mask', '0', '--comp-bias-corr', '0', '--max-seqs', '10000', '-c', '0.8', '-e', '0.001'])
    subprocess.call(['mmseqs', 'filterdb', str(files.res) + '_prof_search',
     str(files.res) + '_' + str(iter_counter) + '_matches_only_db' + '_upd_res','--filter-file', 'target_protID_cluster_file_idx_matches_only_sorted'])
    
    subprocess.call(['mmseqs', 'result2msa', 
    #files.query_db + '_clu' + '_rep' + '_profile',
    files.query_db + '_clu_msa_db_profile',
     files.target_db,
     str(files.res) + '_' + str(iter_counter) + '_matches_only_db' + '_upd_res',
     str(files.res) + '_' + str(iter_counter) + '_matches_only_db' + '_upd_res_msa',
     '--msa-format-mode', '4'])

    subprocess.call(['mmseqs', 'convertmsa', 
     str(files.res) + '_' + str(iter_counter) + '_matches_only_db' + '_upd_res_msa',
     str(files.res) + '_' + str(iter_counter) + '_matches_only_db' + '_upd_res_msa_db'])
    
    subprocess.call(['mmseqs', 'msa2profile', str(files.res) + '_' + str(iter_counter) + '_matches_only_db' + '_upd_res_msa_db', str(files.res) + '_' + str(iter_counter) + '_matches_only_db' + '_upd_res_msa_db_profile'])



def add_new_profiles(clu_indices_for_frac_occ_min_df):    
    # now let's construct the msa for the new proteins from the neighbourhood to add them to the query (concatenate with the msa of the old proteins + matches from above in function make_new_query)
    subprocess.call(['mmseqs', 'cluster', 
    str(files.res) + '_' + str(iter_counter) +'_neigh_only_db',
     str(files.res) + '_' + str(iter_counter) +'_neigh_only_db' + '_clu',
     'tmp', '-c', '0.8', '-e', '0.001'])

    # that's to reject sequence clusters when the fraction of cluster matches in which they occur is smaller than frac_occ_min
    subprocess.call(['mmseqs', 'createtsv', 
    str(files.res) + '_' + str(iter_counter) +'_neigh_only_db',
     str(files.res) + '_' + str(iter_counter) +'_neigh_only_db',
     str(files.res) + '_' + str(iter_counter) +'_neigh_only_db' + '_clu',
     str(files.res) + '_' + str(iter_counter) +'_neigh_only_db' + '_clu.tsv'])
    
    clu_tsv = pd.read_csv(str(files.res) + '_' + str(iter_counter) +'_neigh_only_db' + '_clu.tsv', dtype=None, sep='\t', header = None)
    clu_tsv.columns = ['rep', 'real_prot_id']
    # now let's merge the clustering results df with the df containing counts of spatial clusters per protein
    clu_tsv_spatial_clu_counts = pd.merge(clu_tsv, clu_indices_for_frac_occ_min_df, on='real_prot_id')
    #print(clu_tsv_spatial_clu_counts)
    frac_occ_min = float(args.frac_occ_min)
    print(frac_occ_min)
    # count in how many spatial clusters every profile occurs
    fraction_df = pd.DataFrame()
    print(clu_tsv_spatial_clu_counts)
    #fraction_df = clu_tsv_spatial_clu_counts.groupby(['rep'])['internal_clu_ind'].expanding().count().reset_index(name='spatial_clu_n')
    fraction_df = clu_tsv_spatial_clu_counts.groupby(['rep'])['internal_clu_ind'].transform('count').reset_index(name='spatial_clu_n')
    fraction_df['real_prot_id'] = clu_tsv_spatial_clu_counts['real_prot_id']
    fraction_df['rep'] = clu_tsv_spatial_clu_counts['rep']
    #print(fraction_df)
    total_number_of_spatial_clusters = len(clu_tsv_spatial_clu_counts['internal_clu_ind'].unique())
    #print('total_number_of_spatial_clusters', total_number_of_spatial_clusters)
    fraction_df['fraction'] = fraction_df['spatial_clu_n']/total_number_of_spatial_clusters
    #print(fraction_df)
    prot_filtered_by_fraction = fraction_df.loc[fraction_df['fraction'] > frac_occ_min]['real_prot_id']
    print(prot_filtered_by_fraction)
    neigh_lookup = pd.read_csv(str(files.res) + '_' + str(iter_counter) +'_neigh_only_db.lookup', dtype=None, sep='\t', header = None)
    
    neigh_lookup.columns = ['int_id', 'real_prot_id', 'source_id']
    merged_filter_frac_lookup_df = neigh_lookup.merge(prot_filtered_by_fraction, on='real_prot_id', how='inner')
    print(merged_filter_frac_lookup_df)
    merged_filter_frac_lookup_df['int_id'].to_csv('seq_clu_filter_fraction_idx', index=False, header=False)  

    subprocess.call(['mmseqs', 'createsubdb', 'seq_clu_filter_fraction_idx',
     str(files.res) + '_' + str(iter_counter) +'_neigh_only_db' + '_clu',
     str(files.res) + '_' + str(iter_counter) +'_neigh_only_db' + '_clu' + '_frac_filter'])
    
    subprocess.call(['mmseqs', 'result2msa', 
    str(files.res) + '_' + str(iter_counter) +'_neigh_only_db',
     str(files.res) + '_' + str(iter_counter) +'_neigh_only_db',
     str(files.res) + '_' + str(iter_counter) +'_neigh_only_db' + '_clu' + '_frac_filter',
     str(files.res) + '_' + str(iter_counter) +'_neigh_only_db' + '_clu' + '_frac_filter_msa',
     '--msa-format-mode', '4'])

    subprocess.call(['mmseqs', 'convertmsa', 
     str(files.res) + '_' + str(iter_counter) +'_neigh_only_db' + '_clu_frac_filter_msa',
     str(files.res) + '_' + str(iter_counter) +'_neigh_only_db' + '_clu_frac_filter_msa_db'])
    
    subprocess.call(['mmseqs', 'msa2profile', str(files.res) + '_' + str(iter_counter) +'_neigh_only_db' + '_clu_frac_filter_msa_db', str(files.res) + '_' + str(iter_counter) +'_neigh_only_db' + '_clu_msa_db_profile'])


def make_new_query():
    # the function makes the new query database by concatenating (old updated profiles + new proteins) msas from update_profiles() and add_new_profiles()
    print('making new query db')
    query_db_path = str(files.query_db)
    logging.debug(f"iter_counter {iter_counter}")
    # this idea from the above comment disabled as it causes problem if there are number in the file name
    #if iter_counter > 1:
    #    query_db_path = str(files.query_db)[:(str(files.query_db).find(str(iter_counter-1)))]
    logging.debug(f"query_db_path {query_db_path}")

    subprocess.call(['mmseqs', 'concatdbs', str(files.res) + '_' + str(iter_counter) + '_matches_only_db' + '_upd_res_msa_db_profile', str(files.res) + '_' + str(iter_counter) +'_neigh_only_db' + '_clu_msa_db_profile', query_db_path + str(iter_counter) + 'iter_db_clu_msa_db_profile'])

    subprocess.call(['mmseqs', 'concatdbs', str(files.res) + '_' + str(iter_counter) + '_matches_only_db' + '_upd_res_msa_db_profile_h', str(files.res) + '_' + str(iter_counter) +'_neigh_only_db' + '_clu_msa_db_profile_h', query_db_path + str(iter_counter) + 'iter_db_clu_msa_db_profile_h'])

def set_match_threshold(match_score_gap, query_specific_thresholds):
    # the function calculates the simple heuristic to set profile-specific match score threshold. That is needed to only have positionally ortholodous families
    # match score is the mmseqs2 similarity bit score
    # here to calculate the threshold only for matches that were part of the cluster the filtered results file (only with cluster matches, without single matches) is used
    filtered_search_res_clu_matches_path = str(files.res) + '_' + str(iter_counter) + '_matches_only_db' + '_upd_res'

    subprocess.call(['mmseqs', 'convertalis',
    # files.query_db + '_clu' + '_rep' + '_profile',
    files.query_db + '_clu_msa_db_profile',
     files.target_db, filtered_search_res_clu_matches_path,
      filtered_search_res_clu_matches_path +'.m8'])

    filtered_search_res_only_clu_matches_df = pd.read_csv(str(filtered_search_res_clu_matches_path)+'.m8', dtype={'str':'float'}, sep='\t', header = None)
    filtered_search_res_only_clu_matches_df.columns = ['query_id', '2','3', '4', '5', '6', '7', '8', '9', '10', '11', 'bit_score']
    print('filtered_search_res_clu_matches', filtered_search_res_only_clu_matches_df)
    
    min_scores_per_query_df = filtered_search_res_only_clu_matches_df.loc[filtered_search_res_only_clu_matches_df.groupby('query_id').bit_score.idxmin()][['query_id','bit_score']]
    min_scores_per_query_df['bit_score'] = min_scores_per_query_df['bit_score'] - match_score_gap
    print('min_scores_per_query_df', min_scores_per_query_df)
    print('query_specific_thresholds', query_specific_thresholds)
    dict_match_thresholds_upd = pd.Series(min_scores_per_query_df.bit_score.values,index=min_scores_per_query_df.query_id).to_dict()
    query_specific_thresholds.update(dict_match_thresholds_upd)
    print('query_specific_thresholds', query_specific_thresholds)



def initialize_new_prot_score2(old_query_upd_scores, arr_clu_neigh_prots, arr_matches_in_clu, query_specific_thresholds):
    print('initializing new proteins scores')
    # here i get real indices for all the new proteins, from the neighbourhood and inside of the clusters (non-matches)
    arr_prot_to_add = np.setdiff1d(arr_clu_neigh_prots, arr_matches_in_clu)

    print('arr_prot_to_add', arr_prot_to_add)
    print('length arr_prot_to_add', len(arr_prot_to_add))
    
    # Unlike before, here I just set slightly positive scores for the new proteins from the neighbourhood and non-matches from inside, i dont calculate the concrete scores anymore as i dont do searches anymore
    arr_score_x = np.empty(len(arr_prot_to_add))
    #arr_score_x.fill(0.4)
    arr_score_x.fill(1)
    print('arr_score_x', arr_score_x)
    

    # Here I add just to the scores dict the dict made out of new prots ids and their scores
    #dict_additional_scores = dict(zip(arr_prot_to_add, arr_score_x))
    dict_additional_scores = dict(zip(arr_prot_to_add, arr_score_x))
    old_query_upd_scores.update(dict_additional_scores)
    # also add new threshold (set at first to 0) to the dictionary with query-specific match score thresholds
    thresholds_to_init = np.empty(len(arr_prot_to_add))
    thresholds_to_init.fill(0)
    dict_additional_thresholds = dict(zip(arr_prot_to_add, thresholds_to_init))
    query_specific_thresholds.update(dict_additional_thresholds)

    logging.debug(f"old_query_upd_scores updated with neighbours \n {old_query_upd_scores}")
    return(old_query_upd_scores)


# This part is for "iteration 0", flag which is switched when 
# you have some proportion of single queries not forming clusters
def find_singletons(mapped_res):
    mapped_results = mapped_res.res_map_to_header
    index_list = mapped_res.ind_list
    target_db_lookup = mapped_res.target_db_lookup
    target_db_h = mapped_res.target_db_h

    # For some reasons, _h file is not sorted as lookup by default, so I sort it accordingly
    # here list(set()) is to remove duplicates of genomes/proteins if any
    target_db_h['sort_cat'] = pd.Categorical(target_db_h['ID'], categories=list(set(target_db_lookup.iloc[:, 1].tolist())), ordered=True)
    target_db_h.sort_values('sort_cat', inplace=True)
    target_db_h.reset_index(inplace=True)

    # to fix the problem with the nan coming from reading the table
    logging.debug(f"mapped_results \n {mapped_results}")
    logging.debug(f"index_list \n {index_list}")

    cluster_matches = list()
    matches_ids_list = mapped_results['ID'].tolist()

    score_max_cluster = 1
    
    for i in range(0, len(target_db_lookup.iloc[:, 1])-1):
        if target_db_h["ID"].values[i] in matches_ids_list:
            curr_query_id = mapped_results.loc[mapped_results['ID'] == target_db_h["ID"].values[i], 'query_ID'].iloc[0]
            strand = [int(target_db_h["strand"].values[i])]
            # making list just to fit to the format
            # expected by extract_proteins_neighbourhood
            target_hit = [target_db_h["ID"].values[i]]
            i_0_cluster_start = int(target_db_h["coord1"].values[i])
            i_1_cluster_end = int(target_db_h["coord2"].values[i])
            # here I add twice the left and the right coord (as lists) to make it compatible with the other parts where I have a list for all coordinates of prots
            cluster_matches.append((i_0_cluster_start,
         i_1_cluster_end, score_max_cluster, 
         curr_query_id, target_hit, strand, [i_0_cluster_start], [i_1_cluster_end]))
 
    logging.debug(f"cluster_matches \n {cluster_matches}")
    return cluster_matches


def run_singleton_search():
    # Running search to get if the neighbourhood proteins are enriched 
    # in these "pseudoclusters"
    print('running search for singletons scores')
    subprocess.call(['mmseqs', 'cluster', files.query_db, files.query_db + '_clu',
     'tmp', '--min-seq-id', '0.9'])
    subprocess.call(['mmseqs', 'createsubdb', files.query_db + '_clu', files.query_db,
     files.query_db + '_clu' + '_rep'])
    subprocess.call(['mmseqs', 'createsubdb', files.query_db + '_clu', files.query_db + '_h',
     files.query_db + '_clu' + '_rep' + '_h'])
    subprocess.call(['mmseqs', 'result2profile', files.query_db + '_clu' + '_rep',
     files.query_db, files.query_db + '_clu', files.query_db + '_clu' + '_rep' + '_profile' ])
    search_cov = str(args.search_cov)
    subprocess.call(['mmseqs', 'search', 
    files.query_db + '_clu' + '_rep' + '_profile',
     files.target_db,
     files.res + '_prof_search',
     'tmp', '-a', '--mask', '0', '--comp-bias-corr', '0', '--max-seqs', '10000', '-c', search_cov, '-e', '0.001'])
    subprocess.call(['mmseqs', 'convertalis', files.query_db + '_clu' + '_rep' + '_profile',
     files.target_db, files.res + '_singleton_prof_search',
      files.res + '_singleton_prof_search' +'.m8'])


def initialize_singleton_score(sign_clusters_df, mapped_res):
    # FIX IT!!! setting up the hits score to 1 (as at the begin of the 1st iter)
    # and the neighbourhood prots to 0 currently as it happens practically
    # for neighbourhood scores for other iterations
    print('calculating singleton score')
    #search_result_file = pd.read_csv(files.res + '_singleton_prof_search' +'.m8', dtype={'str':'float'}, sep='\t', header = None)
    query_db_path = str(files.query_db)
    new_query_db_lookup = pd.read_csv(query_db_path+'.lookup', dtype=None, sep='\t', header = None)
    logging.debug(f"{query_db_path+str(iter_counter) + 'iter_db.lookup'}")
    #crass_query_db0iter_db
    mapped_results = mapped_res.res_map_to_header
    old_query_upd_scores = dict()
    
    # without 'to list' it doesnt check if value in in column
    mapped_res_ID_list = mapped_results['ID'].tolist()
    mapped_res_query_ID_list = mapped_results['query_ID'].tolist()
    for new_prot_id in new_query_db_lookup.loc[:,1]:
        #print('new_prot_id', new_prot_id)
        if new_prot_id in mapped_res_ID_list or new_prot_id in mapped_res_query_ID_list:
            old_query_upd_scores[new_prot_id] = 1
        else:
            old_query_upd_scores[new_prot_id] = 0
        logging.debug(f"old_query_upd_scores \n {old_query_upd_scores}")
    return(old_query_upd_scores)


def preprocess_singleton_main(old_query_upd_scores, UpdatedStats):

    s_0 = UpdatedStats['s_0']
    d_strand_flip_penalty = UpdatedStats['d_strand_flip_penalty']


    make_profiles()
    run_search()
    mapped_res = ResultsMapping.map_target_to_coord(query_specific_thresholds)
    mapped_res.res_map_to_header.to_csv('mapped_results_mish', sep = '\t')

    cluster_matches = find_singletons(mapped_res)
    
    cluster_matches_df = pd.DataFrame(cluster_matches)
    cluster_matches_df.to_csv('single_matches_raw', sep = '\t', index = False)
    print('number of singletons', len(cluster_matches_df.index))
    
    sign_clusters_df = cluster_matches_df.copy()
    sign_clusters_df.columns = ["coord1", "coord2", "score",
     "query", "target_prots", "strand", "coordS1", "coordS2"]
    
    logging.debug(f"sign_clusters_df \n {sign_clusters_df}")

    # For singletons, I do not extract proteins from sides, in the function 
    # if_snigleton = '1' the number of prots from 1 side set to 0. So I only extract matches 
    #extract_proteins_cluster_neighborhood(sign_clusters_df, mapped_res)
    #make_new_query()

    #
    
    arr_mmseqs_ind_matches_in_clu, arr_mmseqs_ind_clu_neigh_only, arr_clu_neigh_prots, arr_matches_in_clu, clu_indices_for_frac_occ_min_df = extract_proteins_cluster_neighborhood(sign_clusters_df, mapped_res)
    reassign_non_enriched(old_query_upd_scores, bias, s_0)
    update_profiles()
    add_new_profiles(clu_indices_for_frac_occ_min_df)
    make_new_query()
    set_match_threshold(match_score_gap, query_specific_thresholds)
    old_query_upd_scores = initialize_new_prot_score2(old_query_upd_scores, arr_clu_neigh_prots, arr_matches_in_clu, query_specific_thresholds)

    UpdatedStats['s_0'] = s_0
    UpdatedStats['d_strand_flip_penalty'] = d_strand_flip_penalty

    #old_query_upd_scores = initialize_singleton_score(sign_clusters_df, mapped_res)



#def cluster_clusters(significant_cluster_df_enriched):
def cluster_clusters(significant_cluster_df_enriched):
    print('clustering clusters')
    if args.evalfilteruse == '1':
        path_to = files.res + '_' + str(iter_counter) + '_sign_clusters_df_filter_upd'
    else:
        path_to = files.res + '_' + str(iter_counter) + '_sign_clusters_df_filter_upd'
    path_to_test = path_to
    #path_to_test = '/Users/Sasha/Documents/GitHub/mishpokhe_test/anti_crispr_res_wOld_6_guo1iter_res_2_iter_sign_clusters_enrich_stat_filtered'
    # ast.literal_eval is used to read rows looking like [1,2,3] as python list
    clusters_stat = pd.read_csv(path_to_test, dtype=None, sep='\t',
     converters={'query_prots':ast.literal_eval,'target_prots':ast.literal_eval, 'strand':ast.literal_eval, "coordS1":ast.literal_eval, "coordS2":ast.literal_eval})
    #clusters_stat = significant_cluster_df_enriched
    clusters_queries = clusters_stat['query_prots'].copy()
    queries_set = set()
    for row in clusters_queries:
        for query in row:
            queries_set.add(query)
    try:
        queries_set.remove('')
    except KeyError:
        pass
    queries_list = list(queries_set)
    #print(queries_set)
    #print(len(queries_set))
    # Dont forget that set makes random order, so order in list every run of function varies
    clusters_dict = dict()
    i = 0
    list_presence_lists = list()
    dict_presence_lists = dict()
    for row in clusters_queries:
        query_presence_list = np.zeros(len(queries_set))
        cluster_set = set(row)
        try:
            cluster_set.remove('')
        except KeyError:
            pass
        #print(cluster_set)
        # removing singleton clusters here
        #if len(cluster_set) == 1:
        #    continue
        clusters_dict[i] = row
        list_presence_lists.append(query_presence_list)
        #print(list_presence_lists[i])
        for prot in cluster_set:
            if prot in queries_list:
                ind = queries_list.index(prot)
                list_presence_lists[i][ind] = 1    
        dict_presence_lists[i] = list_presence_lists[i]
        i = i + 1
        #print(queries_list[ind])
        #print(list_presence_lists[i-1])
        #print(dict_presence_lists)
    print(len(list_presence_lists))
    #print(clusters_dict)
    # shape of stacked array = (214, 167) = (queries number, clusters number)
    array_presence_arrays = np.stack(list_presence_lists, axis=0)



    def R_L_density_clustering(dict_presence_lists):
        #print('clustering of clusters start')
        # From https://www.science.org/doi/10.1126/science.1242072
        # clustering from Rodriguez and Laio
        # CHANGE to set cutoff automatically
        # that's to require at least 0.5 of clusters members to be the same, dynamically 
        # changed in the loop
        cutoff_dist = 0.5
        points = dict()
        # making dict of dicts to save distance between data points (spatial clusters)

        @nb.njit(fastmath=True,error_model="numpy")
        #,parallel=True)
        def faster_dist_calc():
            distance_mat_dict = nb.typed.Dict()
            points = nb.typed.Dict()
            for spatial in range(array_presence_arrays.shape[0]):
                print(spatial)
                distance_mat_dict[spatial] = nb.typed.Dict.empty(nb.types.int64, nb.types.float64)
                vec0 = array_presence_arrays[spatial]
                loc_density = 0
                for spatial_compare in range(array_presence_arrays.shape[0]):
                    vec1 = array_presence_arrays[spatial_compare]
                    align_part = np.sum(np.logical_and(vec1,vec0))
                    len_vec1 = np.count_nonzero(vec1 == 1)
                    len_vec0 = np.count_nonzero(vec0 == 1)
                    dist = 1 - (align_part/np.sqrt(len_vec0*len_vec1))
                    distance_mat_dict[spatial][spatial_compare] = dist
                        
                    if dist < cutoff_dist:
                        loc_density = loc_density + 1
                points[spatial] = loc_density
            print(points)
            return points, distance_mat_dict
        points, distance_mat_dict = faster_dist_calc()
        

        # sort points dict by density
        dens_sort_points = dict(sorted(points.items(), key=lambda item: item[1]))
        print(dens_sort_points)
        # get the minimum distance to higher density point
        min_distance_from_higher = dict()
        # dictionary to store the nearest neighbours of higher density for later steps
        closest_higher_dens_neigh = dict()
        # that is a list of spatial cluster ids sorted by density
        list_of_sorted_spatials = list(dens_sort_points.keys())
        print(list_of_sorted_spatials)
        for spatial in list_of_sorted_spatials:
            current_spatial_ind = list_of_sorted_spatials.index(spatial)
            for i in list_of_sorted_spatials[current_spatial_ind+1:len(list_of_sorted_spatials)]:
                if dens_sort_points[i] == dens_sort_points[spatial]:
                    if distance_mat_dict[spatial][i] > cutoff_dist:
                        current_spatial_ind = current_spatial_ind + 1
                    else:
                        break 
                else:
                    break
            print(current_spatial_ind)
            print('current')

            current_comparisons = list_of_sorted_spatials[current_spatial_ind+1:len(list_of_sorted_spatials)]
            # the intermediate dict to find min contains not so many values as many distances
            # are exactly equal and therefore removed from dict, only uniq left
            try:
                min_dist = min({distance_mat_dict[spatial][k] for k in current_comparisons})
                closest_neighbor_high_dens = (list(distance_mat_dict[spatial].keys())[list(distance_mat_dict[spatial].values()).index(min_dist)])
                current_dist_mat_list = [distance_mat_dict[spatial][k] for k in current_comparisons]
                min_dist = min(current_dist_mat_list)
                min_dist_ind = current_dist_mat_list.index(min_dist)
                #closest_neighbor_high_dens = current_comparisons[min_dist_ind]
                ##closest_neighbor_high_dens = (list(distance_mat_dict[spatial].keys())[list(distance_mat_dict[spatial].values())[current_spatial_ind+1:].index(min_dist)])
                if min_dist == 0:
                    pass
                    #q = 2
                    #while min_dist == 0:
                    #    current_comparisons = list_of_sorted_spatials[current_spatial_ind+q:len(list_of_sorted_spatials)]
                    #    min_dist = min({distance_mat_dict[spatial][k] for k in current_comparisons})
                    #    q = q + 1
                #closest_neighbor_high_dens = (list(distance_mat_dict[spatial].keys())[list(distance_mat_dict[spatial].values()).index(min_dist)])
                #closest_neighbor_high_dens = current_comparisons[min_dist_ind]
                #closest_neighbor_high_dens = (list(distance_mat_dict[spatial].keys())[list(distance_mat_dict[spatial].values())[current_spatial_ind+1:].index(min_dist)])

            # that's for the last, highest-dens point, as current_comparisons are empty
            except ValueError:
                min_dist = max(distance_mat_dict[spatial].values())
            #print('current_comparisons', current_comparisons)
            #print('distance_mat_dict[spatial]', distance_mat_dict[spatial])
            #print({distance_mat_dict[spatial][k] for k in current_comparisons})
            print(min_dist)
            min_distance_from_higher[spatial] = min_dist
            #closest_higher_dens_neigh[spatial] = closest_neighbor_high_dens
        print(min_distance_from_higher)
        #for i in closest_higher_dens_neigh.keys():
        #    print(clusters_dict[i], clusters_dict[closest_higher_dens_neigh[i]])


        # measure recommended by paper to find centroids?

        dens_x_dist = dict()
        for k in min_distance_from_higher.keys():
            dens_x_dist[k] = min_distance_from_higher[k]*dens_sort_points[k]

        #import matplotlib.pyplot as plt
        #plt.scatter(dens_sort_points.values(), min_distance_from_higher.values())
        #plt.scatter(dens_x_dist.keys(), dens_x_dist.values())
        #plt.ylim(top=100)
        #plt.ylim(top=500)
        #plt.ylim(top=10)
        #plt.ylim(bottom = 0, top=20)
        #plt.show()
        #print(x)


        all_dist = [distance_mat_dict[k][n]**2 for k in distance_mat_dict.keys() for n in distance_mat_dict[k].keys()]
        #plt.hist(all_dist, bins=range(20))
        #plt.show()

        cluster_centroids = []
        cluster_centroids_ids = []
        # CHECK if I should use dens_x_dist
        dens_x_dist_threshold = 2
        dens_threshold = 400
        dist_threshold = 0.8
        for p in min_distance_from_higher.keys():
            #if min_distance_from_higher[p] > dist_threshold and dens_sort_points[p] < dens_threshold:
            if min_distance_from_higher[p]*dens_sort_points[p] > dens_x_dist_threshold:
            #if min_distance_from_higher[p] > dist_threshold:
                # Ask Johannes if it is correct
                # not letting identical duplicates to become separate cluster centroids
                if clusters_dict[p] not in cluster_centroids:
                    cluster_centroids.append(clusters_dict[p])
                    cluster_centroids_ids.append(p)
        # That happens sometimes, maybe if threshold is not correct
        if not cluster_centroids_ids:
            cluster_centroids_ids = list(min_distance_from_higher.keys())
            cluster_centroids = [clusters_dict[p] for p in min_distance_from_higher.keys()]
        print(cluster_centroids)
        print(len(cluster_centroids))
        print(cluster_centroids_ids)

        final_clusters_ids1 = {k: [] for k in cluster_centroids_ids}
        final_clusters_reals1 = {', '.join(k): [] for k in cluster_centroids}
        print(final_clusters_ids1)
        print(final_clusters_reals1)
        print(dens_sort_points)
        print('assigning points to clusters')
        iterate_over1 = list(dens_sort_points.keys())[::-1]
        for point in iterate_over1:
            print('point', point)
            if point not in final_clusters_ids1:
                distances_to_centroids = {ke: distance_mat_dict[point][ke] for ke in final_clusters_ids1.keys()}
                min_dist = min(distances_to_centroids.values())
                #print('final_clusters_ids.keys()', final_clusters_ids.keys())
                #print('distance_mat_dict[point]', distance_mat_dict[point])
                #print('distances_to_centroids', distances_to_centroids)
                #print(point_dist_mat)
                if min_dist < 1.0:
                    # that is the closest centroid
                    curr_centroid = list(distances_to_centroids.keys())[list(distances_to_centroids.values()).index(min_dist)] 
                    print('curr_centroid', curr_centroid)
                    final_clusters_ids1[curr_centroid].append(point)
                    final_clusters_reals1[', '.join(clusters_dict[curr_centroid])].append(clusters_dict[point])
                else:
                    print('curr centroid is itself')
                    final_clusters_ids1[point] = []
                    final_clusters_reals1[', '.join(clusters_dict[point])] = []
                    cluster_centroids.append(clusters_dict[point])
                    cluster_centroids_ids.append(point)
        print(final_clusters_ids1)
        #print(final_clusters_reals)
        print('here are your clusters')
        raw_clu_of_clu = open(str(files.res) + str(iter_counter) +'_clu_of_clu_all', 'w')
        
        for k in final_clusters_reals1.keys():
            raw_clu_of_clu.write(f'centroid is {str(k)} \n')
            for n in final_clusters_reals1[k]:
                raw_clu_of_clu.write(str(n)+'\n')
            raw_clu_of_clu.write('-------'+'\n')
        print('number of clusters is', len(cluster_centroids))
        #print('intercentroid dist')
        #for point in final_clusters_ids1.keys():
        #    print({ke: distance_mat_dict[point][ke] for ke in final_clusters_ids1.keys()})

        #for c in cluster_centroids:
        #    print(c)
        #    print('**')

        final_clusters_ids = final_clusters_ids1
        final_clusters_reals = final_clusters_reals1 

        dict_presence_lists_incl_clu = copy.deepcopy(dict_presence_lists)
        clu_ids_for_id_file = dict()

        for k in final_clusters_ids.keys():
        #    dict_presence_lists_incl_clu[k] = np.append(dict_presence_lists_incl_clu[k]
            clu_ids_for_id_file[k] = k
            for v in final_clusters_ids[k]:
        #        dict_presence_lists_incl_clu[v] = np.append(dict_presence_lists_incl_clu[v],k)
                clu_ids_for_id_file[v] = k


        print('clu_ids_for_id_file   ')
        print(clu_ids_for_id_file) 
        
        #print(dict_presence_lists_incl_clu)
        # BE CAREFUL! the last values are cluster centroids
        #  Looks like: {....,291: (array([0., 0., 0., 0., 1., ..., 70),
        # 292: (array([0., 0., 0., ..., 0., 292)} 70 and 292 are centroids of the cluster here
        print(final_clusters_ids)
        #print(final_clusters_reals)
        print(len(final_clusters_ids.keys()))
        print('iterate_over1', len(iterate_over1))
        print('dens_sort_points', dens_sort_points)
        return final_clusters_ids, clu_ids_for_id_file
    
    final_clusters_ids, clu_ids_for_id_file = R_L_density_clustering(dict_presence_lists)
    #print(x)



    not_clustered_to_initial_acrs = list()
    clustered_to_initial_acrs = list()

    for l in final_clusters_ids.keys(): 
        #print(l)
        associated = 0
        #print(clusters_dict[l])
        #if 'anti' in ''.join(clusters_dict[l]):
        if clusters_stat['initial_q_or_match'][l] == True:
            associated = 1
        else:
            for i in final_clusters_ids[l]:
                string_of_queries = ''.join(clusters_dict[i])
                #if 'anti' in string_of_queries:
                if clusters_stat['initial_q_or_match'][i] == True:
                    #print(string_of_queries)
                    associated = 1
                    break
        #print('associated is, ', associated)
        if associated == 1:
            clustered_to_initial_acrs.append(l)
            clustered_to_initial_acrs.extend(final_clusters_ids[l])
            #clustered_to_initial_acrs.append('----')
        else:
            not_clustered_to_initial_acrs.append(l)
            not_clustered_to_initial_acrs.extend(final_clusters_ids[l])
    print(clustered_to_initial_acrs)
    positives_filtered = []


    #close_and_associated = list(set(clustered_to_close) & set(clustered_to_initial_acrs))
    # it checks here if list 'close_and_associated' is not empty
    #if iter_counter > 1 and not close_and_associated == False:
    #    significant_clusters_eval_filter_df_clu = clusters_stat.iloc[sorted(close_and_associated)]
    #else:
    #    significant_clusters_eval_filter_df_clu = clusters_stat.iloc[sorted(clustered_to_initial_acrs)]
    ids_column_list = list()
    for i in sorted(list(clu_ids_for_id_file.keys())):
            ids_column_list.append(clu_ids_for_id_file[i])

    clusters_stat["arc_clu_id"] = ids_column_list
    
    significant_clusters_eval_filter_df_clu = clusters_stat.iloc[sorted(clustered_to_initial_acrs)]

    # make a cluster ids column
    
    #clu_arc_ind = open(str(files.res) + str(iter_counter) +'clu_arc_ind', 'w')
    #clu_arc_ind.write('-------'+'\n')
    # do not really understand why sorting needed in the next line. But otherwise it gets errors trying to assess
    # non-existing elements of old_query_scores in the next iter (probably something related to the order?)
    if iter_counter == 2:
        for i in sorted(clustered_to_initial_acrs):
            print(clusters_dict[i])
        #print(x)
    for c in clustered_to_initial_acrs:
        #print('query string', clusters_dict[c])
        print([clusters_stat['target_prots'][c], clusters_stat['coord1'][c], clusters_stat['coord2'][c]])
        positives_filtered.append([clusters_stat['target_prots'][c], clusters_stat['coord1'][c], clusters_stat['coord2'][c]])
    print('filtered positives', len(positives_filtered))

    # --------new clusterting filter
    return significant_clusters_eval_filter_df_clu


def main(old_query_upd_scores, UpdatedStats):

    s_0 = UpdatedStats['s_0']
    d_strand_flip_penalty = UpdatedStats['d_strand_flip_penalty']

    import time
    start_time = time.time()

    #if iter_counter > 1:
    if if_singleton == 1:
        if iter_counter == 0:
            make_profiles()
    else:
        if iter_counter == 1:
            make_profiles()
            
    print("--- %s seconds for make_profiles() ---" % (time.time() - start_time))
    # CHECK if it works with multihitdb (just from command line it worked)
    # CHECK why there are more results with multihitdb (target is converted to profiles??)
    #if iter_counter > 1:
    start_time1 = time.time()
    run_search()
    print("--- %s seconds for run_search() ---" % (time.time() - start_time1))
    # this class is to have order and strand for target proteins
    # FIX best query hit is needed
    
    # real mapping, uncomment? DELETE? as there are duplicates of this command in init

    # MAKE coord and strand integers
    start_time2 = time.time()
    mapped_res = ResultsMapping.map_target_to_coord(query_specific_thresholds)
    mapped_res.res_map_to_header.to_csv('mapped_results_mish', sep = '\t')
    print("--- %s seconds for map_target_to_coord() ---" % (time.time() - start_time2))

    # Is it ok to assign to None?
    start_time3 = time.time()
    use_intermediate = 0
    cluster_matches_fname = str(files.res) + str(iter_counter) + 'cluster_matches'
    if use_intermediate == 1 and iter_counter == 1:
        f=open(cluster_matches_fname,"r")
        lst=f.read()
        f.close()
        cluster_matches=eval(lst)
    else:
        cluster_matches = find_clusters(mapped_res, old_query_upd_scores, d_strand_flip_penalty, s_0, enrichment_bias)
        f=open(cluster_matches_fname,"w")
        f.write(str(cluster_matches))
        f.close()
    print("--- %s seconds for find_clusters() ---" % (time.time() - start_time3))
    
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    
    logging.debug(f'checking global {len(old_query_upd_scores)} {d_strand_flip_penalty}, {s_0}')
    cluster_matches_df = pd.DataFrame(cluster_matches)
    cluster_matches_df.to_csv(str(files.res) + str(iter_counter) + 'cluster_matches_raw', sep = '\t', index = False)
    print('number of clusters', len(cluster_matches_df.index))

    # That is to keep intermediate cluster matches files
    cluster_matches_fname = str(files.res) + str(iter_counter) + 'cluster_matches'
    f=open(cluster_matches_fname,"w")
    f.write(str(cluster_matches))
    f.close()

    print(cluster_matches_df)
    # CHECK if it is optimal to divide scores update to few functions

    # REMOVE duplicated function calling?
    # REMOVE these functions?
    #update_scores_for_cluster_matches(cluster_matches)
    #significant_cluster_df_enriched, s_0 = update_scores_for_cluster_matches(cluster_matches, mapped_res)
    #print(significant_cluster_df_enriched)
    
    # UNCOMMENT
    start_time4 = time.time()
    stat_lambda, stat_K = calculate_karlin_stat(cluster_matches, mapped_res, s_0, bias)
    print("--- %s seconds for calculate_karlin_stat() ---" % (time.time() - start_time4))
    start_time5 = time.time()
    significant_cluster_df, significant_clusters_eval_filter_df = calculate_e_value(stat_lambda, stat_K, cluster_matches, mapped_res)
    print("--- %s seconds for calculate_e_value() ---" % (time.time() - start_time5))
    start_time6 = time.time()
    significant_cluster_df_enriched, s_0, old_query_upd_scores, L, l = update_scores_for_cluster_matches(significant_clusters_eval_filter_df, mapped_res, bias)
    print("--- %s seconds for update_scores_for_cluster_matches() ---" % (time.time() - start_time6))
    sign_clusters_df = significant_cluster_df_enriched
    
    # Checking if user wants to use architecture clustering-filtering
    arc_filter = int(args.arc_filter)
    if arc_filter == 1:
        start_time7 = time.time()
        significant_clusters_eval_filter_df_clu = cluster_clusters(significant_cluster_df_enriched)
        print("--- %s seconds for cluster_clusters() ---" % (time.time() - start_time7))
    #if iter_counter == 2:
    #    print(x)

    # CHANGE notation for significant clusters??

    d_strand_flip_penalty = set_strand_flip_penalty(cluster_matches, mapped_res)

    # CHECK do I need this step??
    #sign_clusters_df = set_strand_flip_penalty(cluster_matches)

    if arc_filter == 1: 
        # Is it okay to use e-val filtering here?
        sign_clusters_df = significant_clusters_eval_filter_df_clu
        ###sign_clusters_df = significant_cluster_df_enriched

        path_clu_filter = files.res + '_' + str(iter_counter) + '_iter_sign_clusters_enrich_stat_filtered_clu_filter'
        significant_clusters_eval_filter_df_clu.to_csv(path_clu_filter, sep = '\t', index = False)


    arr_mmseqs_ind_matches_in_clu, arr_mmseqs_ind_clu_neigh_only, arr_clu_neigh_prots, arr_matches_in_clu, clu_indices_for_frac_occ_min_df = extract_proteins_cluster_neighborhood(sign_clusters_df, mapped_res)
    reassign_non_enriched(old_query_upd_scores, bias, s_0)
    update_profiles()
    add_new_profiles(clu_indices_for_frac_occ_min_df)
    make_new_query()
    set_match_threshold(match_score_gap, query_specific_thresholds)
    old_query_upd_scores = initialize_new_prot_score2(old_query_upd_scores, arr_clu_neigh_prots, arr_matches_in_clu, query_specific_thresholds)

    UpdatedStats['s_0'] = s_0
    UpdatedStats['d_strand_flip_penalty'] = d_strand_flip_penalty


    #generate_mmseqs_ffindex(sign_clusters_df)


if __name__ == "__main__":

    print("starting")
    arg_parser()
    iterations = int(args.iter)
    if_singleton = int(args.singleton)

    # To write log file from scratch, otherwise would add to the file content
    if os.path.exists(f"{args.res}_mishpokhe.log"):
        os.remove(f"{args.res}_mishpokhe.log")
    logging.basicConfig(filename=f"{args.res}_mishpokhe.log", format='%(name)s - %(levelname)s -%(funcName)s() - %(message)s', level=logging.DEBUG)
    
    files = FilePath.get_path()
    print(files.query_db)
    print(files.target_db)
    print(files.res)

    # get queries ids to fill the scores dict with 1
    sep = r"'\t'"
    awk_print = "'{print $2}'"
    cmd = f'awk -F {sep} {awk_print} {files.query_db}.lookup'
    #print(cmd)
    q_and_matches = list()
    # WARNING: popen is deprecated, think about it
    queries = os.popen(cmd).read().splitlines()
    old_query_upd_scores = dict()
    query_specific_thresholds = dict()
    for q in queries:
        #print(q)
        old_query_upd_scores[q] = 1
        query_specific_thresholds[q] = 0
        #print(old_query_upd_scores)

    # that is to write all the INITIAL queries to this list which allow later filter in the clustering of clusters. The filter would keep those clusters  which are clustered together with those having match to the initial query/another match to the initial query
    query_db_lookup = pd.read_csv(str(files.query_db)+str(".lookup"), dtype=None, sep='\t', header = None)
    q_and_matches = query_db_lookup.iloc[:, 1].tolist()

    # FIX to be the right order of functions (should be after run_search())

    # UNCOMMENT, real mapping!
    #mapped_res = ResultsMapping.map_target_to_coord()

    # Make changeable
    bias = 0 #4
    enrichment_bias = 0 #4
    match_threshold = 0
    match_score_gap = 10

    #Keeps and updates d_strand_flip_penalty and s_0
    UpdatedStats = dict()
    UpdatedStats['d_strand_flip_penalty'] = None
    UpdatedStats['s_0'] = -1 #-4.4 #-0.4 - bias

    # For 0th iteration with query containing singletons
    print(if_singleton)
    if if_singleton == 1:
        print('doing singletons')
        iter_counter = 0
        # that is to not mix 0 and 1 iter res
        saved_files_res = files.res
        files.res = str(files.res) + '0iter'
        preprocess_singleton_main(old_query_upd_scores, UpdatedStats)
    
        files.query_db = str(files.query_db) + str(iter_counter) + 'iter_db'
        files.res = saved_files_res

    iter_counter = 1
    while iterations > 0:
        print('iter_counter:',iter_counter)
        logging.debug(f"iter_counter, files.query_db: {iter_counter, files.query_db}")
        if iter_counter == 1:
            files.query_db = files.query_db
        if iter_counter == 2:
            files.query_db = str(files.query_db) + str(iter_counter-1) + 'iter_db'
            files.res = str(files.res) + str(iter_counter-1) + 'iter_res'
        if iter_counter > 2:
            # The enhancement below is disabled for now as might cause problem if the filepath has other numbers in the name
            #query_db_path = str(files.query_db)[:str(files.query_db).find(str(iter_counter-2))]
            query_db_path = str(files.query_db)
            files.query_db = query_db_path + str(iter_counter-1) + 'iter_db'
            res_path = str(files.res)[:str(files.res).rfind(str(iter_counter-2))]
            files.res = res_path + str(iter_counter-1) + 'iter_res'
        logging.debug(f"files.query_db in main: {files.query_db}")

        main(old_query_upd_scores, UpdatedStats)
        iterations -= 1
        iter_counter += 1
    
    # Here it is to merge the results from all iteration to one file (I have to take care of the competitive things as well)


