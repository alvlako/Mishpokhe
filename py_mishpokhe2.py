from cmath import inf
import ctypes as ct
from ctypes import *
import numpy as np
import pandas as pd

import itertools
import os
import re
import subprocess
import sys
from sys import stdout

# version 2 accordingly to the written in the proposal

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
                query_db = input('Enter query_db path: ')
                target_db = input('Enter target_db path: ')
                res = input('Enter results file path: ')
                return self(query_db,target_db, res)
            except:
                print('Invalid input!')
                continue


# add parameters adjustment option? error raising?
# deleting previous tmp and db files
def make_profiles():
    print('making query mmseqs profiles')
    subprocess.call(['mmseqs', 'cluster', files.query_db, files.query_db + '_clu',
     'tmp', '--min-seq-id', '0.9'])
    subprocess.call(['mmseqs', 'createsubdb', files.query_db + '_clu', files.query_db,
     files.query_db + '_clu' + '_rep'])
    subprocess.call(['mmseqs', 'createsubdb', files.query_db + '_clu', files.query_db + '_h',
     files.query_db + '_clu' + '_rep' + '_h'])
    subprocess.call(['mmseqs', 'result2profile', files.query_db + '_clu' + '_rep',
     files.query_db, files.query_db + '_clu', files.query_db + '_clu' + '_rep' + '_profile' ])
    pass


# add options?
# CHECK if its right to search query profiles against target
# ADD besthit
def run_search():
    print('running mmseqs profile search')
    # extracting best hit per match with --max-accept 1
    # ASK Johannes if it is ok to use this without sensitivity iters (manual)
    # Ask Johannes if it is correct to search TARGET against Query
    # FIX, --max-accept 1 does not work for now
    subprocess.call(['mmseqs', 'search', 
    files.query_db + '_clu' + '_rep' + '_profile',
     files.target_db,
     files.res + '_prof_search',
     'tmp', '-a'])
    #subprocess.call(['mmseqs', 'search', files.target_db,
    #files.query_db + '_clu' + '_rep' + '_profile',
    # files.res + '_prof_search',
    # 'tmp', '-a'])
    #, '--max-accept', '1'])
    # convert to convenient format
    subprocess.call(['mmseqs', 'convertalis', files.query_db + '_clu' + '_rep' + '_profile',
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
     def map_target_to_coord(self):
        print('mapping results')
        #search_result_file = pd.read_csv(str(files.res + '_prof_search'), dtype={'int':'float'}, sep='\t')
        #search_result_file = pd.read_csv('/Users/Sasha/Documents/GitHub/mishpokhe_test/py2_multihit_res', dtype={'int':'float'}, sep='\t')
        #search_result_file = pd.read_csv('/Users/Sasha/Documents/GitHub/mishpokhe_test/py_norm_res_prof_search.m8', dtype={'str':'float'}, sep='\t', header = None)
        search_result_file = pd.read_csv(files.res + '_prof_search' +'.m8', dtype={'str':'float'}, sep='\t', header = None)
        print(search_result_file)
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
        print(search_result_file)
        #print(search_result_file.groupby('target_id')["eval"].transform('min'))
        search_result_file = search_result_file.loc[search_result_file.groupby('target_id')['eval'].idxmin()].reset_index(drop=True)

        # REMOVE unique? it's unique from the prev step
        real_id_list = pd.Series(pd.unique(search_result_file.iloc[:, 1]))
        print("real ids list", real_id_list)
        # map by 1st (0) column with real ids from search res
        # print(target_db_h.loc[target_db_h.iloc[:, 0].astype(str) == 'MT006214.1_1'])
        tmp_res_map_to_header = target_db_h.loc[target_db_h.iloc[:, 0].astype(str).isin(real_id_list)]


        # NOTe - isin makes sorting, so the next lines to get to the original lines order
        # that is to link proper query ids to the other info
        tmp_res_map_to_header['sort_cat'] = pd.Categorical(tmp_res_map_to_header['ID'], categories=real_id_list, ordered=True)
        tmp_res_map_to_header.sort_values('sort_cat', inplace=True)
        tmp_res_map_to_header.reset_index(inplace=True)
        tmp_res_map_to_header['query_ID'] = search_result_file.iloc[:, 0]

        # sort it back again to iterate via by order in the next step
        # these ids are strings, need to be numeric for correct sorting
        
        def atoi(text):
            return int(text) if text.isdigit() else text

        def natural_keys(text):
            return [ atoi(c) for c in re.split(r'(\d+)', text) ]
        
        print(' ')
        print(' ')
        print(' ')
        print('start soritng')
        print(' ')

        human_ids = tmp_res_map_to_header['ID'].tolist()
        human_ids.sort(key=natural_keys)
        #tmp_res_map_to_header['ID'] = human_ids.astype('category')
        tmp_res_map_to_header['ID_cat'] = pd.Categorical(tmp_res_map_to_header['ID'], categories=human_ids, ordered=True)
        #tmp_res_map_to_header.sort_values('ID_cat')
        
        tmp_res_map_to_header.sort_values('ID_cat', inplace=True)
        tmp_res_map_to_header.reset_index(inplace=True, drop=True)
        res_map_to_header = tmp_res_map_to_header.copy()

        #res_map_to_header = tmp_res_map_to_header.sort_values(by=['ID'])

        ##res_map_to_header['ind'] = ind_list.values

        # sorting target db lookup to iterate correctly in find_clusters()
        # DOUBLE check?
        human_ids2 = target_db_lookup['ID'].tolist()
        human_ids2.sort(key=natural_keys)
        target_db_lookup['ID_cat'] = pd.Categorical(target_db_lookup['ID'], categories=human_ids2, ordered=True)
        target_db_lookup.sort_values('ID_cat', inplace=True)
        target_db_lookup.reset_index(inplace=True, drop=True)
        target_db_lookup = target_db_lookup.drop(['ID_cat'], axis=1)


        # get target proteins real ids

        # CLEAN
        ind_list = real_id_list

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
def find_clusters(mapped_res, old_query_upd_scores, d_strand_flip_penalty, s_0):

    mapped_results = mapped_res.res_map_to_header
    results = mapped_res.search_result_file
    index_list = mapped_res.ind_list
    target_db_lookup = mapped_res.target_db_lookup
    target_db_h = mapped_res.target_db_h

    # For some reasons, _h file is not sorted as lookup by default, so I sort it accordingly
    target_db_h['sort_cat'] = pd.Categorical(target_db_h['ID'], categories=target_db_lookup.iloc[:, 1].tolist(), ordered=True)
    target_db_h.sort_values('sort_cat', inplace=True)
    target_db_h.reset_index(inplace=True)

    # to fix the problem with the nan coming from reading the table
    results = results[results.iloc[:, 0].notna()]
    print(results)
    print('mapped results')
    print(mapped_results)
    print("index list", index_list)

    #Algorithm 1 - iterate through target prot
    print(results.iloc[:, 0].size)

    cluster_matches = list()
    # CHECK if score max cluster set up correct
    score_max_cluster = 0
    # OPTIMIZE how to set the strand penalty
    if d_strand_flip_penalty is None:
        d_strand_flip_penalty = 1
    else:
        d_strand_flip_penalty = -1*d_strand_flip_penalty
    
    # CHECK if it was set up correctly
    score_min_cluster = 0
    score_i_minus_1_cluster = 0

    if s_0 is None:
        #s_0 = -0.01
        s_0 = -1

    # Shows if the current protein has a query match
    is_match = 1

    # CHECK if correct (esp if cluster does not start from the 1st gene)
    i_0_cluster_start = int(target_db_h["coord1"].values[0])
    pd.set_option('display.max_columns', None)
    #print(mapped_results)
    print("cluster start", i_0_cluster_start)
    i_1_cluster_end = int(target_db_h["coord2"].values[0])
    # CHECK if correct
    init_strand = int(target_db_h["strand"].values[0])
    print('strand is')
    print(init_strand)

    query_genes_ids = []
    target_genes_ids = []
    prots_strands = []


    matches_ids_list = mapped_results['ID'].tolist()

    # CHECK it now ignores the last line, is it a constant feature to
    # have empty line at the ened of the results file???
    # FIX genes ids retrieval
    print('target lookup')
    print(target_db_lookup.iloc[:, 1])
    for i in range(0, len(target_db_lookup.iloc[:, 1])-1):
        print('NEW cycle starts')
        print(i)
        print(target_db_lookup.iloc[:, 1].values[i])
        diff_genomes_penalty = 0

        
        # CHANGE this score (to 0 and 1 for 1st iter?)
        #score_x_i = float(results.iloc[i,10])
        # CHECK if correct
        if target_db_h["ID"].values[i] in matches_ids_list:
            score_x_i = 1
            is_match = 1
            curr_query_id = mapped_results.loc[mapped_results['ID'] == target_db_h["ID"].values[i], 'query_ID'].iloc[0]
            if old_query_upd_scores is not None:
                score_x_i = old_query_upd_scores[curr_query_id]
        else:
            score_x_i = s_0
            is_match = 0
            # is it a good idea?
            curr_query_id = ''
        # FIX temporary solution to make score_x_i to overweight other scores to get significant clusters
        #score_x_i = -np.log(score_x_i)
        # CHECK in evalue section how to initialize this score
        # CHECK if -1 (and other potential values) properly read
        strand = int(target_db_h["strand"].values[i])
        print("strand is", strand)

        # gap changed to use gene number, not the coordinate diff
        #gap = abs(int(mapped_results["coord1"].values[i])- int((mapped_results["coord2"].values[i-1]))) - 1
        # CHECK if i need 1 or 2 column after split (for now it is for id looks like NC_015295.1_78)

        # I need no gap penalty as I have s0
        #print(target_db_h["ID"].values[i].split("_"))
        #gap = abs(int(target_db_h["ID"].values[i].split("_")[1])- int((target_db_h["ID"].values[i-1].split("_")[1])))
        #print('gap =', gap)
        
        # CHECK if that's strand of the previous gene or next
        # FIX with OR statement (re.findall('\+|\-', test))
        # FIX something wrong the first str = init gives "different strand"
        # FIX use another, calculated strand flip penalty (d) in the next iterations
        if strand == init_strand:
            f_strand_flip = 0
            print('same strand')
        else:
            print('different strand')
            f_strand_flip = 1

        # to check whether proteins are from the same genome
        # THINK if it should be done better
        # also it relies on having "." in prot id
        # CHECK if I should compare with the prev prot
        # ADD something because now it depends a lot on format, indices would be different here
        # if ids are NC_111.1_1 or something else like MT333.1_2
        # Here for the format of  NC_111.1_1
        #print(target_db_h["ID"].values[i].split("_")[0]+target_db_h["ID"].values[i].split("_")[1])
        #print(target_db_h["ID"].values[i+1].split("_")[0]+target_db_h["ID"].values[i+1].split("_")[1])
        # Here for the format of MT333.1_2
        #print(target_db_h["ID"].values[i].split("_")[0])
        #print(target_db_h["ID"].values[i+1].split("_")[0])
        # CHANGE here if the dataset big or small
        # Here for the format of MT333.1_2
        #if target_db_h["ID"].values[i].split("_")[0]!= target_db_h["ID"].values[i+1].split("_")[0]:
        # Here for the format of  NC_111.1_1
        if target_db_h["ID"].values[i].split("_")[0]+target_db_h["ID"].values[i].split("_")[1] != target_db_h["ID"].values[i+1].split("_")[0]+target_db_h["ID"].values[i+1].split("_")[1]:
            print('different genomes!!!')
            # is it a good idea?
            diff_genomes_penalty = 10000
            #continue

        #print('passed')
        # updating previous gene strand (current gene = previous for next for loop iter)
        init_strand= strand
        
        print('scores')
        print('prev cluster score', score_i_minus_1_cluster)
        print('strand flip', f_strand_flip*d_strand_flip_penalty)
        print('current score', score_x_i)
        print('  ')
        print(score_i_minus_1_cluster - f_strand_flip*d_strand_flip_penalty + score_x_i)
        print(max(0, score_x_i))


        if (score_i_minus_1_cluster - f_strand_flip*d_strand_flip_penalty - diff_genomes_penalty + score_x_i) > max(0, score_x_i):
            score_i_cluster = score_i_minus_1_cluster - f_strand_flip*d_strand_flip_penalty  - diff_genomes_penalty + score_x_i
            print('proceed', score_i_cluster, score_max_cluster)
            print('score i-1', score_i_minus_1_cluster)
            
            if score_i_cluster > score_max_cluster:
                print('first')
                score_max_cluster = score_i_cluster
                # CHECK if correct, CHANGE to look better
                #i_1_cluster_end = str(mapped_results.iloc[i,1])
                i_1_cluster_end = int(target_db_h["coord2"].values[i])
                # CHECK changing of score_i_minus_1_cluster
                # CHECK why did add this? it is not in the Johannes pseudocode
                #score_i_minus_1_cluster = score_i_cluster
                print('upd prev cluster score', score_i_minus_1_cluster)

                #last_target_gene = mapped_results["ID"].values[i]
                # ADD best match??
                print('first 1st append')
                query_genes_ids.append(curr_query_id)
                target_genes_ids.append(target_db_h["ID"].values[i])
                prots_strands.append(int(target_db_h["strand"].values[i]))

        else:
            print('second')
            score_i_cluster = score_i_minus_1_cluster - f_strand_flip*d_strand_flip_penalty  - diff_genomes_penalty + score_x_i
            print('second_proceed', score_i_cluster, score_max_cluster)
            print('score i-1', score_i_minus_1_cluster)
            score_i_cluster = score_x_i

            # CHECK if correct, CHANGE to get right coord
            # CHECK, maybe I have to take i-1 coord? otherwise
            # it seems like I take the left coord of 1st non-cluster
            # prot - ASK Johannes
            # THINK if the next line here or below
            #i_0_cluster_start = int(mapped_results["coord1"].values[i])
            
            print(score_max_cluster, score_min_cluster)
            if score_max_cluster > score_min_cluster:
                print('second 1st append')
                cluster_matches.append((i_0_cluster_start,
                i_1_cluster_end, score_max_cluster,
                query_genes_ids, target_genes_ids, prots_strands))
                score_max_cluster = 0
            # THINK if it is okay to be here or above
            # for some reasons if here it gives proper result
            i_0_cluster_start = int(target_db_h["coord1"].values[i])

            #first_target_gene = mapped_results["ID"].values[i]
            query_genes_ids = []
            # ADD best match??
            query_genes_ids.append(curr_query_id)
            target_genes_ids = []
            target_genes_ids.append(target_db_h["ID"].values[i])
            prots_strands = []
            prots_strands.append(int(target_db_h["strand"].values[i]))
            #query_genes_ids.append(mapped_results["query_ID"].values[i])
            #target_genes_ids.append(mapped_results["ID"].values[i])
        # CHECK if correct, not as in latex
        score_i_minus_1_cluster = score_i_cluster
        print('max and min scores', score_max_cluster, score_min_cluster)
        print('cluster coord', i_0_cluster_start, i_1_cluster_end)
        #print('cluster matches', cluster_matches)

    if score_max_cluster > score_min_cluster:
        print('second 2nd append')
        cluster_matches.append((i_0_cluster_start,
         i_1_cluster_end, score_max_cluster, 
         query_genes_ids, target_genes_ids, prots_strands))
    # add more to cluster matches table? prot id?
    # ADD return of the changed ResultsMapping object? (with added scores?)
    # FIX to be faster or remove
    print(cluster_matches)
    return cluster_matches


# it re-defines the scores which I have got in the first iteration 
# ADD significant clusters determination?
def update_scores_for_cluster_matches(cluster_matches, mapped_res):
    # CHANGE for really significant clusters, not as it is?
    significant_clusters = cluster_matches 
    # ADD query id to mapped results
    mapped_results = mapped_res.res_map_to_header
    results = mapped_res.search_result_file

    target_db_lookup = mapped_res.target_db_lookup
    target_db_h = mapped_res.target_db_h

    # CHECK if these are right columns
    # CHECK if the query and target assignment is correct
    K = len(significant_clusters)
    L = len(target_db_lookup.index)

    print(L)
    
    bias = 0
    sign_clusters_df = pd.DataFrame(significant_clusters)
    sign_clusters_df.columns = ["coord1", "coord2", "score",
    "query_prots", "target_prots", "strand"]

    print(sign_clusters_df)
    sign_clusters_df['new_score_enrich'] = 0
    sign_clusters_df['list_new_scoreS_enrich'] = 0

    # MAKE faster?
    sign_clusters_df['queries_string'] = [','.join(map(str, l)) for l in sign_clusters_df['query_prots']]

    # Should cluster prots be done better?
    cluster_prots = pd.DataFrame()
    cluster_prots['query_id'] = sign_clusters_df['query_prots'].explode()
    cluster_prots['target_id'] = sign_clusters_df['target_prots'].explode()

    l = len(cluster_prots)

    # CHANGE later, THINK
    bias = 0

    old_query_upd_scores = dict()

    print("K, L, l", K, L, l)
    
    print('sign_clusters_df')
    print(sign_clusters_df)
    print('cluster_prots')
    print(cluster_prots)
    print('mapped results')
    print(mapped_results)
    # CHECK if it is correct to iterate through cluster matches?
    # So I am calculating enrichment score for cluster matches and for cluster prots with no match
    # CHECK if leaving only unique entries is correct
    # pseudocounts are added with parameter aplha_pseudocount
    aplha_pseudocount = pow(10,np.log10(1/L)-1)
    x_number_of_queries = len(mapped_results['query_ID'].unique())
    print('x_number_of_queries', x_number_of_queries)
    # Do I need unique?
    for query_id in cluster_prots['query_id'].unique():
        print(query_id)
        # to avoid empty queries corresponding to prots with no match
        if query_id == '':
            continue
        M_x = mapped_results['query_ID'][mapped_results['query_ID'] == query_id].shape[0]
        m_x = cluster_prots[cluster_prots['query_id'] == query_id]['query_id'].count()
        # in this case M_x and m_x are equal as there were no single hits. 
        # CHECK with other data where would be hits not only in clusters
        print("M_x, m_x", M_x, m_x)
        # using pseudocounts for m_x/l proportion to avoid zeros in log
        cluster_prot_proportion = np.divide(m_x, l)
        cluster_prot_proportion = np.divide((m_x+aplha_pseudocount), (l + x_number_of_queries*aplha_pseudocount))
        score_x = np.log(np.divide(cluster_prot_proportion, np.divide(M_x, L))) - bias
        print("updated score for q", query_id, "is", score_x)

        old_query_upd_scores[query_id] = score_x

        # MAKE faster?
        # adding score of the query prot to get summarized score for the cluster
        # CHECK if correct
        sign_clusters_df.loc[sign_clusters_df['queries_string'].str.contains(str(query_id)),
        'new_score_enrich'] = sign_clusters_df.loc[sign_clusters_df['queries_string'].str.contains(str(query_id)),
        'new_score_enrich'] + score_x
        sign_clusters_df.loc[sign_clusters_df['queries_string'].str.contains(str(query_id)),
        'list_new_scoreS_enrich'] = sign_clusters_df.loc[sign_clusters_df['queries_string'].str.contains(str(query_id)),
        'list_new_scoreS_enrich'].astype(str) + ',' + str(score_x)
        #print(sign_clusters_df['list_new_scoreS_enrich'])
        # FIX figure out how to set b and threshold for prot to be enriched in cluster
        # ADD file to store scores corresponding to profiles

    # calculating the scores for cluster prots with no match
    # THINK if it is ok to iterate again and also the order of scores adding to the
    # list of scores is not competely correct
    matches_ids_list = mapped_results['ID'].tolist()

    # Counting all target cluster prots having query matches
    m_x_sum = len(cluster_prots['query_id'])
    # Counting all target proteins with matches to some queries
    M_x_sum = len(mapped_results['ID'])

    for target_id in cluster_prots['target_id']:
        print('the cycle is running')
        if target_id not in matches_ids_list:
            print('not in list')
            print('no match', target_id)
            s_0 = np.divide(np.divide((l-m_x_sum),l), np.divide((L-M_x_sum),L)) - bias
            print('updates s_0 for', target_id, 'is = ', s_0)
            print(sign_clusters_df[pd.DataFrame(sign_clusters_df['target_prots'].tolist()).isin(target_id.split()).any(1).values])
            tmp_df = sign_clusters_df[pd.DataFrame(sign_clusters_df['target_prots'].tolist()).isin(target_id.split()).any(1).values]['new_score_enrich'] + s_0
            array_of_bool = pd.DataFrame(sign_clusters_df['target_prots'].tolist()).isin(target_id.split()).any(1).values
            index_of_row = int(np.where(array_of_bool == True)[0])
            sign_clusters_df.loc[sign_clusters_df.index[index_of_row],'new_score_enrich'] = sign_clusters_df.loc[sign_clusters_df.index[index_of_row],'new_score_enrich'] + s_0
            # FIX to make the cell string 
            #sign_clusters_df.loc[sign_clusters_df.index[index_of_row],'list_new_scoreS_enrich'] = sign_clusters_df.loc[sign_clusters_df.index[index_of_row],'list_new_scoreS_enrich'] + ',' + str(s_0)
            print(sign_clusters_df)
            print(sign_clusters_df['list_new_scoreS_enrich'].tolist())
        else:
            print('is in list')
            s_0 = -1
    if len(cluster_prots['target_id']) == 0:
        print('0 target in clusters')
        s_0 = -1
    print(s_0)
            

    print(sign_clusters_df)
    sign_clusters_df.to_csv('sign_clusters_df', sep = '\t')
    return(sign_clusters_df, s_0, old_query_upd_scores, L, l)


# CHANGE to not duplicate so much the func above
def calculate_karlin_stat(cluster_matches, mapped_res):
    print('using c code for Karlin-Altschul statistics')

    subprocess.call(['cc', '-fPIC', '-shared', '-o', 'karlin_c.so', 'karlin_c.c'])
    # !! CHANGE the path later
    so_file = "./karlin_c.so"
    my_functions = CDLL(so_file)

    # ASK Johannes if that's ok that I calculated enrich scores for all search results,
    # not only the clusters

    mapped_results = mapped_res.res_map_to_header

    print('calculating prob for each match')
    print('mapped results')
    print(mapped_results)

    all_enrich_scores = []

    significant_clusters = cluster_matches 

    # CHECK if these are right columns
    # CHECK if the query and target assignment is correct
    K = len(significant_clusters)
    
    target_db_lookup = mapped_res.target_db_lookup
    L = len(target_db_lookup.index)
    
    bias = 0
    sign_clusters_df = pd.DataFrame(significant_clusters)
    sign_clusters_df.columns = ["coord1", "coord2", "score",
    "query_prots", "target_prots", "strand"]

    cluster_prots = pd.DataFrame()
    cluster_prots['query_id'] = sign_clusters_df['query_prots'].explode()
    cluster_prots['target_id'] = sign_clusters_df['target_prots'].explode()

    l = len(cluster_prots)
    print('len(cluster_prots)', l)
    aplha_pseudocount = pow(10,np.log10(1/L)-1)
    
    # CHECK if it is right that in case of all res probs I should iterate through target
    # NOT through query as in example above
    n_target = len(mapped_results['ID'].unique())
    for target_id in mapped_results['ID'].unique():
        print(target_id)
        M_x = mapped_results['ID'][mapped_results['ID'] == target_id].shape[0]
        m_x = cluster_prots[cluster_prots['target_id'] == target_id]['target_id'].count()
        # in this case M_x and m_x are equal as there were no single hits. 
        # CHECK with other data where would be hits not only in clusters
        print("M_x, m_x", M_x, m_x)
        # using pseudocounts for m_x/l proportion to avoid zeros in log
        # ASK Johannes if pseudocount is correct
        cluster_prot_proportion = np.divide(m_x, l)
        cluster_prot_proportion = np.divide((m_x+aplha_pseudocount), (l + n_target*aplha_pseudocount))
        score_x = np.log(np.divide(cluster_prot_proportion, np.divide(M_x, L))) - bias
        print("updated score for t", target_id, "is", score_x)
        all_enrich_scores.append(score_x)
    
    print('scores for all res', all_enrich_scores)



    enrich_scores = np.array(all_enrich_scores)
    print('len', len(enrich_scores))

    #if 0 not in enrich_scores:
    #    no_0 = True
    #    enrich_scores = np.append(enrich_scores, 0)
    #else:
    #    no_0 = False
    #print('no_0', no_0)

    # REMOVE later, just for current test


    print(enrich_scores)
    #print(x)
    # ASK Johannes if I done the scores correctly
     # # expand unique scores in order to have larger range of integer, wider range to calc e-value
    # # multipliying param for log-scores
    #mult_param = 6
    mult_param = int(input('enter multiplying value '))
    scores_bins = np.histogram_bin_edges(enrich_scores, bins='auto')
    print('BINS', scores_bins)
    print('bins diff', scores_bins[1]-scores_bins[0])
    #print(x)
    # # THINK!

    # ASK Johannes, REMOVE later, that's just to fix problem with log when m_x = 0
    # ASK Johannes if it is ok to use binning, no I think, it is NOT ok
    unique_scores = mult_param*np.unique(enrich_scores)
    unique_scores = np.sort(unique_scores)
    print('multiplied, sorted', unique_scores)
    #print(x)

    unique_scores = np.round(unique_scores, decimals=0)
    unique_scores, score_counts = np.unique(unique_scores, return_counts=True)
    unique_scores = unique_scores.astype(int)
    print(unique_scores)
    print('uniq enrich scores', np.unique(enrich_scores))
    print(np.unique(mult_param*enrich_scores))
    print(np.round(unique_scores, decimals=0))
    print('len of orig scores', len(np.unique(enrich_scores)))
    print("len of scores", len(unique_scores))

 

   

    print(score_counts)
    print('len', len(enrich_scores))
    
    score_prob = score_counts / len(enrich_scores)
    print(score_prob)

    # figure out why round works weird so -3.5 > 4, -2.5 > 2
    scores_table = np.column_stack((np.round(unique_scores, decimals=0), score_prob))
    print(scores_table)
    #print(x)
    # MAKE faster!
    # the loop to insert missed integers (should be one by one) to scores and 0 to correspond. probs
    prev_int_val = int(scores_table[0][0])
    curr_row_ind = 0
    for i in scores_table:
        print('start inserting')
        print('int_val', prev_int_val)
        print('i', i)
        # CHECK what negative value in times_insert what do (seems like does nothing correctly)
        times_insert = int(i[0]) - prev_int_val - 1
        print('times_insert', times_insert)
        
        for j in range(0, times_insert):
            insert_index = np.where(np.all(scores_table==i,axis=1))[0][0]
            print('insert_index', insert_index)
            int_to_insert = prev_int_val + j + 1
            scores_table = np.insert(scores_table, insert_index, np.array((int_to_insert, 0)), 0)  
            j = j + 1
        prev_int_val = int(i[0])
        curr_row_ind = curr_row_ind + 1
        print(scores_table)

    # REMOVE!! tmp! for little dataset testing to have 0
    #scores_table = np.insert(scores_table, 1, np.array((0, 0.33)), 0)  
    #scores_table = np.insert(scores_table, 2, np.array((1, 0.17)), 0)  
    #print(scores_table)
    
    unique_scores = scores_table[:,0]
    score_prob = scores_table[:,1]

    # ASK Johannes what to do if all the scores are negative?
    unique_scores
    index_of_0 = np.where(unique_scores == 0)[0][0]
    
    
    print('final uniq and prob', unique_scores, score_prob)
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
    print('min and max', min_score, max_score)

    
    print("index_of_0", index_of_0)
    print('prob of 0', arr[index_of_0:])
    print('arr',arr)
    arr = np.array(arr)
    print(type(arr))
    print('element type', type(arr[0]))
    
    pointer_to_middle = arr[index_of_0:].ctypes.data_as(ct.POINTER(ct.c_double))
    pointer_to_middle = arr[index_of_0:].ctypes.data_as(ct.POINTER(ct.c_double))
    print(pointer_to_middle.contents)

    
    my_functions.BlastKarlinLambdaNR.argtypes = [POINTER(c_double),c_int32, c_int32]
    my_functions.BlastKarlinLambdaNR.restype = c_double
    BlastKarlin_lambda = my_functions.BlastKarlinLambdaNR(pointer_to_middle, min_score, max_score)
    print("lambda", BlastKarlin_lambda)
    BlastKarlin_lambda = ct.c_double(BlastKarlin_lambda)
    
    my_functions.BlastKarlinLtoH.restype = c_double
    BlastKarlin_H = my_functions.BlastKarlinLtoH(pointer_to_middle, min_score, max_score,
     BlastKarlin_lambda)
    print('H', BlastKarlin_H)
    BlastKarlin_H = ct.c_double(BlastKarlin_H)
    

    my_functions.BlastKarlinLHtoK.restype = c_double
    BlastKarlin_K = my_functions.BlastKarlinLHtoK(pointer_to_middle, min_score, max_score,
     BlastKarlin_lambda, BlastKarlin_H)
    print('K', BlastKarlin_K)
    return(BlastKarlin_lambda, BlastKarlin_K)



# ASK where strand flip penalty goes? next iter?
def set_strand_flip_penalty(cluster_matches, mapped_res):
    target_db_h = mapped_res.target_db_h
    mapped_results = mapped_res.res_map_to_header
    #print(target_db_h)
    # this just delete the consecutive duplicates, therefore I can cound 
    # how many switches between strands are in the file (1 1 -1 -1 1 -> 1 -1 1)
    F = len([key for key, _group in itertools.groupby(target_db_h["strand"])]) - 1
    print("F = ", F)
    # CHANGE to not duplicate update_scores
    significant_clusters = cluster_matches
    sign_clusters_df = pd.DataFrame(significant_clusters)
    sign_clusters_df.columns = ["coord1", "coord2", "score",
    "query_prots", "target_prots", "strand"]
    print(sign_clusters_df)
    f = 0
    # this should count number of flips in each cluster
    # think how to speed up this step
    for cluster_strands in sign_clusters_df["strand"]:
        n_flips_cluster = len([key for key, _group in itertools.groupby(cluster_strands)]) - 1
        f = f + n_flips_cluster
    print("f = ", f)
    # CHANGE to not duplicate scores updating
    K = len(significant_clusters)
    cluster_prots = pd.DataFrame()
    cluster_prots['query_id'] = sign_clusters_df['query_prots'].explode()
    l = len(cluster_prots)
    print("K, l = ", K, l)
    target_db_lookup = mapped_res.target_db_lookup
    L = len(target_db_lookup.index)

    # ASK Johannes if pseudocount is correct
    aplha_pseudocount = pow(10,np.log10(1/L)-1)
    # CHECK if it is right that in case of all res probs I should iterate through target
    # NOT through query as in example above
    n_target = len(mapped_results['ID'].unique())

    #if f == 0:
    #    f = F/100000
    cluster_flips_proportion = np.divide(f, (l - K))
    cluster_flips_proportion = np.divide(f + aplha_pseudocount, (l + n_target*aplha_pseudocount - K))
    d = np.log(np.divide(cluster_flips_proportion, np.divide(F, (L - 1))))
    print(d)
    # is it a good place to return?
    # CHECK if these are really significant clusters
    return(d)


# ASK Johannes if in e-value calculations the L should actually be the L, not like L*L
def calculate_e_value(stat_lambda, stat_K, significant_cluster_df_enriched, mapped_res):
    # FIX to be variable taken from number of prots in target
    # CHANGE it to be linked and the same with the used in calculate_karlin_stat
    mult_param = 6
    mult_param = int(input('enter multiplying value '))
    print('calculating e-value')
    target_db_lookup = mapped_res.target_db_lookup
    L = len(target_db_lookup.index)
    print(significant_cluster_df_enriched)
    stat_lambda = stat_lambda.value
    print(stat_lambda, stat_K)
    significant_cluster_df_enriched["e-value"] = 0
    # PUT in correct place
    pd.options.display.float_format = '{:,.4f}'.format
    for r in range(0, len(significant_cluster_df_enriched.index)):
        print(r)
        enrich_score = significant_cluster_df_enriched["new_score_enrich"][r]
        # REMOVE if not needed
        enrich_scores_list = significant_cluster_df_enriched["list_new_scoreS_enrich"][r].split(',')
        # CHANGE to delete 0 happened with creation of pandas column
        del(enrich_scores_list[0])
        # converting enrich score as for karlin stats calculation
        conv_enrich_score = int(np.round(enrich_score*mult_param, decimals = 0))
        print('conv_enrich_score', conv_enrich_score)
        # is it okay if that's a sum?
        evalue = stat_K*L*np.exp(stat_lambda*(-1)*conv_enrich_score)
        significant_cluster_df_enriched.iat[r, significant_cluster_df_enriched.columns.get_loc("e-value")] = float(evalue)
        print('eval =', evalue)
    # CHANGE the name? 
    # they are not filtered by e-val yet
    cluster_path = files.res + '_' + str(iter_counter) + '_iter_sign_clusters_enrich_stat'
    significant_cluster_df_enriched.to_csv(cluster_path, sep = '\t', index = False)
    # Check if e-val = 0.01 is good filtering
    eval_filter = 1
    use_eval_filter = 1
    if use_eval_filter == 1:
        significant_clusters_eval_filter_df = significant_cluster_df_enriched.loc[(significant_cluster_df_enriched['e-value'] <eval_filter) & (significant_cluster_df_enriched['e-value'] > 0)]
        cluster_path2 = files.res + '_' + str(iter_counter) + '_iter_sign_clusters_enrich_stat_filtered'
        significant_clusters_eval_filter_df.to_csv(cluster_path2, sep = '\t', index = False)
    else:
        print('e-value filter disabled')
        significant_clusters_eval_filter_df = significant_cluster_df_enriched.copy()
    print('significant_cluster_df_enriched', significant_cluster_df_enriched)
    print(significant_clusters_eval_filter_df)
    return significant_cluster_df_enriched, significant_clusters_eval_filter_df



# THINK about writing already to a file in subprocess Popen
def extract_proteins_cluster_neighborhood(sign_clusters_df):
    print('updating query profile started (extracting prots_')
    # CHECK if only the significant clusters are used
    # getting the target prots matched to query
    print('significant (?) clusters table')
    print(sign_clusters_df)
    
    # calling initial fasta for target and query
    # CHANGE it to not ask user again
    #query_fasta = input('Enter query_sequences (fasta) path: ')
    target_fasta = input('Enter target_sequence (fasta) path: ')

    # adding matched target prots to the query set to make the profiles again
    # THINK if it is a good way to add all matched target prots to all queries and run clustering again

    # extracting the left and the right edges proteins for clusters
    print(sign_clusters_df['target_prots'])
    # extracting proteins matched, within and in +3 left and right neighbourhood of the clusters
    ##target_clusters_within = open("target_clusters_within", "w")
    neighbourhood_path = files.res + '_' + str(iter_counter) + '_iter_target_clusters_neighbourhood'
    
    target_clusters_neighbourhood = open(neighbourhood_path, "w")
    target_clusters_matches = open("target_clusters_matches", "w")
    for target_prot_cluster in sign_clusters_df['target_prots']:
        # MAKE faster
        # that is to save cluster prots ids and grep with them to find cluster matches sequences and 
        # exclude them from new within cluster prots
        # WHY do I need this file?
        target_protID_cluster_file = open("target_protID_cluster_file", "w")
        print(*target_prot_cluster)
        # MAKE faster
        print('grepping matches')
        for protID in target_prot_cluster:
            print(protID)
            # whitespace added to grep e.g. 1_1 but not 1_12, 1_13 etc
            target_protID_cluster_file.write(protID + ' ' + '\n')
            matches_seqs = subprocess.Popen(['grep', '-A1', protID + ' ',
             target_fasta], stdout=subprocess.PIPE)
            seqs, error = matches_seqs.communicate()
            target_clusters_matches.write(seqs.decode("utf-8"))
        # REMOVE this file if you wont grep file in file
        target_protID_cluster_file.close()
        # extracting sequences for matches in cluster and writing to the file
        # file - file grep - would not use for now
        #matches_seqs = subprocess.Popen(['grep', '-A1','-f', 'target_protID_cluster_file',
        # target_fasta], stdout=subprocess.PIPE)
        #seqs, error = matches_seqs.communicate()
        #target_clusters_matches.write(seqs.decode("utf-8"))

        target_prot_left = target_prot_cluster[0]
        target_prot_right = target_prot_cluster[-1]
        print('left and right proteins per cycle')
        print(target_prot_left, target_prot_right)
        # extracting header + seq lines between left and right edge proteins in cluster
        # all prots within cluster extracted, but the rightest dont have seq line after this step, only header
        within = subprocess.Popen(['sed', '-n', '/%s /,/%s /p' %(target_prot_left, target_prot_right) + ' ', target_fasta], stdout=subprocess.PIPE)
        
        # FIX FIX FIX
        #remov_pat = '-- "^--$"'
        #print(remov_pat)
        #p_matches_in_cluster = subprocess.Popen(['grep', '-A1', '-F', '-f', 'target_protID_cluster_file', target_fasta, ], stdout=subprocess.PIPE, encoding='utf-8')
        #p2_matches_in_cluster = subprocess.Popen(['grep', '-v', '--'], stdin=p_matches_in_cluster, stdout=subprocess.PIPE, encoding='utf-8')
        #p_matches_in_cluster.stdout.close()
        #print(("the commandline is {}".format(p_matches_in_cluster.args)))
        #print('these are matches in cluster')
        #matches_in_cluster, err_cl = p2_matches_in_cluster.communicate()
        #print(err_cl)
        #print(matches_in_cluster)



        target_clusters_within, error = within.communicate()
        print(target_clusters_within.decode('utf-8'))

        p1_within = subprocess.Popen(['grep', '-A1', target_prot_right + ' ', target_fasta], stdout=subprocess.PIPE)
        # dont communicate if you dont want the pipe to be closed
        p2 = subprocess.Popen(['grep', '-v', target_prot_right], stdin=p1_within.stdout, stdout=subprocess.PIPE)
        p1_within.stdout.close()
        right_protein_seq, err = p2.communicate()
        print('rest',right_protein_seq.decode('utf-8'))
        #target_clusters_within = within
        # extracting 3 prots before and after + seq line for the last prot of the cluster
        # pipes are to get rid of rightest and leftest prots headers
        #subprocess.call(['grep', '-A7', target_prot_right + ' ', target_fasta], stdout=target_clusters_neighbourhood)
        # filtering only concrete genome in db, CHECK if correct for all formats
        # Does it limit the formats option? So I somehow rely on the fact that I have a protein ID 
        # of the format “MW067000.1_2”, and genome ID is of the format “MW067000.1”.
        #genome_id = target_prot_right.split('_')[0]
        group_sep = ''.join(['--', ' ', '"', '^--$','"'])
        # Changed genome id to format of prot = NC_055116.1_25
        # MAKE more general format!
        # For the format of NC_055116.1_25
        genome_id = target_prot_right.split('_')[0] + '_' + target_prot_right.split('_')[1]
        # For the format of YN055116.1_25
        #genome_id = target_prot_right.split('_')[0]

        # extracting not 3 but 6 prots from each side
        p1_left = subprocess.Popen(['grep', '-A1', genome_id, target_fasta], stdout=subprocess.PIPE)
        p4_left = subprocess.Popen(['grep ' + '-v '+ group_sep], stdin=p1_left.stdout, stdout=subprocess.PIPE, shell = True)
        p2_left = subprocess.Popen(['grep', '-B12', target_prot_left + ' '], stdin=p4_left.stdout, stdout=subprocess.PIPE)
        # dont communicate if you dont want the pipe to be closed
        p3_left = subprocess.Popen(['grep', '-v', target_prot_left], stdin=p2_left.stdout, stdout=subprocess.PIPE)
        p1_left.stdout.close()
        p2_left.stdout.close()
        p4_left.stdout.close()
        output_left,err_left = p3_left.communicate()
        print("left prots")
        print(output_left.decode("utf-8"))
        # counting how many proteins I got (in case if it is the end or the start of the file)
        left_prots_extracted_n = output_left.decode("utf-8").count('>')
        print(left_prots_extracted_n)
        # CHANGE later for variable how many prots to extract? 
        # -1 is added to use it later in tail command
        left_prots_remained_to_extract_n_lines = str((6 - left_prots_extracted_n)*2*(-1))
        print("remained", left_prots_remained_to_extract_n_lines)
        # extract remaining number of proteins from another end of the file if any
        if int(left_prots_remained_to_extract_n_lines)*(-1) > 0:
            print('extracting from the end of file')
            print('genome_id', genome_id)
            p1_from_end = subprocess.Popen(['grep', '-A1', genome_id, target_fasta], stdout=subprocess.PIPE)
            p4_from_end = subprocess.Popen(['grep ' + '-v '+ group_sep], stdin=p1_from_end.stdout, stdout=subprocess.PIPE, shell = True)
            p2_from_end = subprocess.Popen(['tail', left_prots_remained_to_extract_n_lines], stdin=p4_from_end.stdout, stdout=subprocess.PIPE)
            output_file_end,err_file_end = p2_from_end.communicate()
            print('end of the file')
            print(left_prots_remained_to_extract_n_lines)
            print(output_file_end.decode("utf-8"))
            p1_from_end.stdout.close()
            target_clusters_neighbourhood.write(output_file_end.decode("utf-8"))


        p1_right = subprocess.Popen(['grep', '-A1', genome_id, target_fasta], stdout=subprocess.PIPE)
        p2 = subprocess.Popen(['grep ' + '-v '+ group_sep], stdin=p1_right.stdout, stdout=subprocess.PIPE, shell = True)
        p3 = subprocess.Popen(['grep', '-A13', target_prot_right + ' '], stdin=p2.stdout, stdout=subprocess.PIPE)
        p4 = subprocess.Popen(['grep', '-v', target_prot_right], stdin=p3.stdout, stdout=subprocess.PIPE)
        p5 = subprocess.Popen(['tail', '-12'], stdin=p4.stdout, stdout=subprocess.PIPE)
        output_right,err_right = p5.communicate()
        p1_right.stdout.close()
        p4.stdout.close()
        p3.stdout.close()
        #print(output_right.decode("utf-8"))

        right_prots_extracted_n = output_right.decode("utf-8").count('>')
        print(right_prots_extracted_n)
        right_prots_remained_to_extract_n_lines = str((6 - right_prots_extracted_n)*2*(-1))
        # extract remaining number of proteins from another end of the file if any
        if right_prots_remained_to_extract_n_lines != '0':
            p1_from_begin = subprocess.Popen(['grep', genome_id, target_fasta], stdout=subprocess.PIPE)
            p2_from_begin = subprocess.Popen(['head',  right_prots_remained_to_extract_n_lines, target_fasta], stdin=p1_from_begin.stdout, stdout=subprocess.PIPE)
            output_file_begin,err_file_begin = p2_from_begin.communicate()
            p1_from_begin.stdout.close()
            target_clusters_neighbourhood.write(output_file_begin.decode("utf-8"))


        
        target_clusters_neighbourhood.write(output_left.decode("utf-8"))
        target_clusters_neighbourhood.write(output_right.decode("utf-8"))
        target_clusters_neighbourhood.write(target_clusters_within.decode("utf-8"))
        target_clusters_neighbourhood.write(right_protein_seq.decode("utf-8"))
        
        # Add target matches to their neighbourhood
        #subprocess.call(['cat', 'target_clusters_matches', '>>', 'target_clusters_neighbourhood'])
        print('neighbourhood_path', neighbourhood_path)
        command_cat = 'cat ' + 'target_clusters_matches ' + '>> ' + neighbourhood_path
        print('command_cat', command_cat)
        os.system(command_cat)

        # FIX!
        #target_clusters_matches.write(matches_in_cluster)
        
    ##target_clusters_within.close()
    target_clusters_neighbourhood.close()
    # FIX!
    target_clusters_matches.close()
    
    # CHECK if you can make faster
    # CHECK if you need to delete duplicates from target prots
    

    #making a new query db
    # THINK if I want not to renew the file but make a separate file for each iteration
    ##subprocess.call(['mmseqs', 'createdb', 'query_fasta_copy_iter', 'query_fasta_copy_iter_db'])

    pass


# function to initialize score for new protein profiles from neighbourhood
# should potentially be called in update_query_profiles_add_proteins function?
# ASK if I need to search again and do everything like for the matches to get the enrich scores?
# ASK if I can just do it in the next iter with all the other profiles
# ASK if now I calculate the enrichment scores correctly
# ASK Johannes what order should it be done, if I initialize scores for proteins or profiles
def make_new_query():
    print('making new query set')
    print('it is iteration ', iter_counter, ' do -1')
    query_fasta = input('Enter query_sequences (fasta) path: ')
    # to have names like query_prot_db2iter_db, not like query_prot_db2iter_db1iter_db
    query_db_path = str(files.query_db)
    print('iter_counter', iter_counter)
    if iter_counter > 1:
        query_db_path = str(files.query_db)[:(str(files.query_db).find(str(iter_counter-1)))]
    print('query_db_path', query_db_path)
    out_file = query_db_path + str(iter_counter) + 'iter.fasta'
    print('out_file', out_file)
    neighbourhood_path = files.res + '_' + str(iter_counter) + '_iter_target_clusters_neighbourhood'
    command_cat = 'cat ' + query_fasta + ' ' + neighbourhood_path + ' >' + ' ' + out_file
    print('command_cat', command_cat)
    os.system(command_cat)
    #with open(out_file, "w") as query_file:
    #    subprocess.Popen(['cat', query_fasta, 'target_clusters_neighbourhood'], stdout = query_file)
    subprocess.call(['mmseqs', 'createdb', query_db_path + str(iter_counter) + 'iter.fasta', query_db_path + str(iter_counter) + 'iter_db'])

def initialize_new_prot_score2(sign_clusters_df, old_query_upd_scores, L, l, mapped_res):
    query_db_path = str(files.query_db)
    if iter_counter > 1:
        query_db_path = str(files.query_db)[:str(files.query_db).find(str(iter_counter-1))]
    new_query_db_lookup = pd.read_csv(query_db_path+str(iter_counter) + 'iter_db.lookup', dtype=None, sep='\t', header = None)
    mapped_results = mapped_res.res_map_to_header
    cluster_prots = pd.DataFrame()
    cluster_prots['target_id'] = sign_clusters_df['target_prots'].explode()

    bias = 0

    aplha_pseudocount = pow(10,np.log10(1/L)-1)
    # Is it correct to use for pseudocounts?
    # CHANGE to not include query prots to counting??
    x_number_of_new_prots = new_query_db_lookup.loc[:,1].size
    print('x_number_of_new_prots', x_number_of_new_prots)
    for new_prot_id in new_query_db_lookup.loc[:,1]:

        M_x = mapped_results['ID'][mapped_results['ID'] == new_prot_id].shape[0]
        m_x = cluster_prots[cluster_prots['target_id'] == new_prot_id]['target_id'].count()
        print(new_prot_id, M_x, m_x)
        cluster_prot_proportion = np.divide((m_x+aplha_pseudocount), (l + x_number_of_new_prots*aplha_pseudocount))
        score_x = np.log(np.divide(cluster_prot_proportion, np.divide(M_x, L))) - bias
        print("updated score for new prot", new_prot_id, "is", score_x)
        
        if new_prot_id not in old_query_upd_scores.keys():
            # FIX? it happens for proteins from neighbourhood, not from the clusters
            if score_x == inf:
                score_x = 0
            old_query_upd_scores[new_prot_id] = score_x
    
    print(old_query_upd_scores)
    return(old_query_upd_scores)




# This part is for "iteration 0", flag which is switched when 
# you have some proportion of single queries not forming clusters
def find_singletons(mapped_res):
    mapped_results = mapped_res.res_map_to_header
    results = mapped_res.search_result_file
    index_list = mapped_res.ind_list
    target_db_lookup = mapped_res.target_db_lookup
    target_db_h = mapped_res.target_db_h

    # For some reasons, _h file is not sorted as lookup by default, so I sort it accordingly
    target_db_h['sort_cat'] = pd.Categorical(target_db_h['ID'], categories=target_db_lookup.iloc[:, 1].tolist(), ordered=True)
    target_db_h.sort_values('sort_cat', inplace=True)
    target_db_h.reset_index(inplace=True)

    # to fix the problem with the nan coming from reading the table
    results = results[results.iloc[:, 0].notna()]
    print(results)
    print('mapped results')
    print(mapped_results)
    print("index list", index_list)

    cluster_matches = list()
    matches_ids_list = mapped_results['ID'].tolist()

    score_max_cluster = 1
    
    for i in range(0, len(target_db_lookup.iloc[:, 1])-1):
        if target_db_h["ID"].values[i] in matches_ids_list:
            curr_query_id = mapped_results.loc[mapped_results['ID'] == target_db_h["ID"].values[i], 'query_ID'].iloc[0]
            strand = int(target_db_h["strand"].values[i])
            target_hit = target_db_h["ID"].values[i]
            i_0_cluster_start = int(target_db_h["coord1"].values[i])
            i_1_cluster_end = int(target_db_h["coord2"].values[i])
            cluster_matches.append((i_0_cluster_start,
         i_1_cluster_end, score_max_cluster, 
         curr_query_id, target_hit, strand))
 
    print(cluster_matches)
    return cluster_matches


def initialize_singleton_score(sign_clusters_df, mapped_res):
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
    subprocess.call(['mmseqs', 'search', 
    files.query_db + '_clu' + '_rep' + '_profile',
     files.target_db,
     files.res + '_prof_search',
     'tmp', '-a'])
    subprocess.call(['mmseqs', 'convertalis', files.query_db + '_clu' + '_rep' + '_profile',
     files.target_db, files.res + '_prof_search',
      files.res + '_prof_search' +'.m8'])

    # Calculating scores


def preprocess_singleton_main():
    make_profiles()
    run_search()
    mapped_res = ResultsMapping.map_target_to_coord()
    mapped_res.res_map_to_header.to_csv('mapped_results_mish', sep = '\t')

    cluster_matches = find_singletons(mapped_res)
    
    cluster_matches_df = pd.DataFrame(cluster_matches)
    cluster_matches_df.to_csv('single_matches_raw', sep = '\t', index = False)
    print('number of singletons', len(cluster_matches_df.index))

    print(cluster_matches_df)
    
    sign_clusters_df = cluster_matches_df.copy()
    sign_clusters_df.columns = ["coord1", "coord2", "score",
     "query", "target", "strand"]

    extract_proteins_cluster_neighborhood(sign_clusters_df)
    make_new_query()

    files.query_db = str(files.query_db) + str(iter_counter) + 'iter_db'
    old_query_upd_scores = initialize_singleton_score(sign_clusters_df, mapped_res)
    print(x)
    pass


def main(old_query_upd_scores, d_strand_flip_penalty, s_0):
    make_profiles()

    # CHECK if it works with multihitdb (just from command line it worked)
    # CHECK why there are more results with multihitdb (target is converted to profiles??)
    run_search()

    # this class is to have order and strand for target proteins
    # FIX best query hit is needed
    
    # real mapping, uncomment? DELETE? as there are duplicates of this command in init

    # MAKE coord and strand integers

    mapped_res = ResultsMapping.map_target_to_coord()
    mapped_res.res_map_to_header.to_csv('mapped_results_mish', sep = '\t')

    # REMOVE, tmp!
    #target_db_lookup = mapped_res.target_db_lookup
    #pd.set_option('display.max_rows', None)
    #print(target_db_lookup)
    #print(x)
    #for i in range(0, len(target_db_lookup.iloc[:, 1])-1):
    #    print(target_db_lookup.iloc[:, 1].values[i])

    # REMOVE, tmp!
    #mapped_res.res_map_to_header = pd.read_csv('mapped_results_mish', sep = '\t', header = 0)
    #print(mapped_res.res_map_to_header)

    # Is it ok to assign to None?

    cluster_matches = find_clusters(mapped_res, old_query_upd_scores, d_strand_flip_penalty, s_0)

    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    
    cluster_matches_df = pd.DataFrame(cluster_matches)
    cluster_matches_df.to_csv('cluster_matches_raw', sep = '\t', index = False)
    print('number of clusters', len(cluster_matches_df.index))

    print(cluster_matches_df)
    # CHECK if it is optimal to divide scores update to few functions

    # REMOVE duplicated function calling?
    # REMOVE these functions?
    #update_scores_for_cluster_matches(cluster_matches)
    #significant_cluster_df_enriched, s_0 = update_scores_for_cluster_matches(cluster_matches, mapped_res)
    #print(significant_cluster_df_enriched)
    
    # UNCOMMENT
    significant_cluster_df_enriched, s_0, old_query_upd_scores, L, l = update_scores_for_cluster_matches(cluster_matches, mapped_res)
    stat_lambda, stat_K = calculate_karlin_stat(cluster_matches, mapped_res)
    sign_clusters_df = significant_cluster_df_enriched
    significant_cluster_df_enriched, significant_clusters_eval_filter_df = calculate_e_value(stat_lambda, stat_K, significant_cluster_df_enriched, mapped_res)

    # CHANGE notation for significant clusters??

    d_strand_flip_penalty = set_strand_flip_penalty(cluster_matches, mapped_res)

    # CHECK do I need this step??
    #sign_clusters_df = set_strand_flip_penalty(cluster_matches)

    # Is it okay to use e-val filtering here?
    sign_clusters_df = significant_clusters_eval_filter_df

    extract_proteins_cluster_neighborhood(sign_clusters_df)
    make_new_query()
    old_query_upd_scores = initialize_new_prot_score2(sign_clusters_df, old_query_upd_scores, L, l, mapped_res)


    #generate_mmseqs_ffindex(sign_clusters_df)


if __name__ == "__main__":

    print("starting")

    print('type number of iterations')
    #iterations = int(input())
    iterations = int(input())

    print('do you have singleton queries? 1 for yes, 0 for no')
    if_singleton = int(input())
    # ADD check if file exists!!!
    files = FilePath.get_path()
    print(files.query_db)
    print(files.target_db)
    print(files.res)

    # FIX to be the right order of functions (should be after run_search())

    # UNCOMMENT, real mapping!
    #mapped_res = ResultsMapping.map_target_to_coord()

    # For 0th iteration with query containing singletons
    if if_singleton == 1:
        iter_counter = 0
        preprocess_singleton_main()
        

    
    iter_counter = 1
    while iterations > 0:
        print(iter_counter)
        print(files.query_db)
        if iter_counter == 1:
            old_query_upd_scores = None
            s_0 = None
            d_strand_flip_penalty = None
            files.query_db = files.query_db
        if iter_counter == 2:
            files.query_db = str(files.query_db) + str(iter_counter-1) + 'iter_db'
            files.res = str(files.res) + str(iter_counter-1) + 'iter_res'
        if iter_counter > 2:
            query_db_path = str(files.query_db)[:str(files.query_db).find(str(iter_counter-2))]
            files.query_db = query_db_path + str(iter_counter-1) + 'iter_db'
            res_path = str(files.res)[:str(files.res).find(str(iter_counter-2))]
            files.res = res_path + str(iter_counter-1) + 'iter_res'
        print('files.query_db in main', files.query_db)


        #print('here is main running')

        #query_db_path = str(files.query_db)
        #print('iter_counter', iter_counter)
        #if iter_counter > 1:
        #    query_db_path = str(files.query_db)[:(str(files.query_db).find(str(iter_counter-1)))]
        #print('query_db_path', query_db_path)
        #out_file = query_db_path + str(iter_counter) + 'iter.fasta'
        #print('out_file', out_file)
        #print('@1', query_db_path + str(iter_counter) + 'iter.fasta')
        #print('@2', query_db_path + str(iter_counter) + 'iter_db')


        main(old_query_upd_scores, d_strand_flip_penalty, s_0)
        iterations -= 1
        iter_counter += 1


