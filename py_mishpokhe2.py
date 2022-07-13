import numpy as np
import pandas as pd

import os
import re
import subprocess

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
def run_search():
    print('running mmseqs profile search')
    subprocess.call(['mmseqs', 'search', 
    files.query_db + '_clu' + '_rep' + '_profile',
     files.target_db,
     files.res + '_prof_search',
     'tmp', '-a'])
    # convert to convenient format
    subprocess.call(['mmseqs', 'convertalis', files.query_db + '_clu' + '_rep' + '_profile',
     files.target_db, files.res + '_prof_search',
      files.res + '_prof_search' +'.m8'])


# The n of prots in db matches the n of gff CDS CDS!!! entries
# FIX best query hit is needed
# CHECK if the mapping correct 
# CHECK if it quaranteed that indices in results file are same order as in lookup (also in order) 
# AND in _h and the order is always the same
# FIX !!!!!! some real ids are read as NaN
# currently it is search used the normal db for target
# CHECK if it is target proteins, not query
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
        search_result_file = pd.read_csv('/Users/Sasha/Documents/GitHub/mishpokhe_test/py_norm_res_prof_search.m8', dtype={'str':'float'}, sep='\t')
        print(search_result_file)
        # Do i need lookup?
        target_db_lookup = np.genfromtxt(str(files.target_db)+str(".lookup"), dtype = None, delimiter="\t", encoding=None)
        target_db_h = pd.read_csv(str(files.target_db)+str("_h"), sep=' # ',header=None)
        #print(search_result_file.iloc[:, 0])
        #print(target_db_h.iloc[:, 0])
        ##ind_list = search_result_file.iloc[:, 0].dropna().astype(int)
        #print(ind_list)
        # get target proteins real ids
        real_id_list = search_result_file.iloc[:, 1]
        print(real_id_list)
        # map by 1st (0) column with real ids from search res
        #print(target_db_h.loc[target_db_h.iloc[:, 0].astype(str) == 'MT006214.1_1'])
        res_map_to_header = target_db_h.loc[target_db_h.iloc[:, 0].astype(str).isin(real_id_list)]
        ##res_map_to_header['ind'] = ind_list.values
        # get target proteins real ids

        print(res_map_to_header)
        ##return self(search_result_file, target_db_lookup, target_db_h, res_map_to_header, ind_list)
        pass



# ?? is it correct that L is set not of all the target proteins
# but only those which were matched and got to results???
# FIX !!!!!! some real ids are read as NaN
# FIX ! to process not the initial results but the best hit results
# CHECK whether you need to change the scores
# DELETE header from pandas dataframes from ResultMapping class object
# Set the scores not to e-values column but to 0
def find_clusters():

    # FIX order the results by target prot ID
    mapped_results = mapped_res.res_map_to_header
    results = mapped_res.search_result_file
    index_list = mapped_res.ind_list
    # to fix the problem with the nan coming from reading the table
    results = results[results.iloc[:, 0].notna()]
    print(results)
    print('mapped results')
    print(mapped_results)
    print(index_list)

    #Algorithm 1 - iterate through target prot
    print(results.iloc[:, 0].size)

    cluster_matches = list()
    # CHECK if score max cluster set up correct
    score_max_cluster = 0
    # OPTIMIZE how to set the strand penalty
    d_strand_flip_penalty = 0.5
    # CHECK if it was set up correctly
    score_min_cluster = 0
    score_i_minus_1_cluster = 0.1

    # CHECK, ASK Johannes
    gap_penalty = 100

    # CHECK if correct (esp if cluster does not start from the 1st gene)
    i_0_cluster_start = str(mapped_results.iloc[0,1])[:str(mapped_results.iloc[0,1]).find(''.join(re.findall('\+|\-', str(mapped_results.iloc[0,1]))))]
    i_1_cluster_end = 0
    # CHECK if correct
    strand = str(mapped_results.iloc[0,1])
    print('strand is')
    print(strand)
    if re.findall('\+', strand)==[]:
        init_strand= '-'
    if re.findall('\-', strand)==[]:
        init_strand= '+'
    print(init_strand)

    # CHECK it now ignores the last line, is it a constant feature to
    # have empty line at the ened of the results file???
    for i in range(0, len(results.iloc[:, 0])):
        print(i)
        #print(results.iloc[[i]])
        score_x_i = float(results.iloc[i,3])
        # FIX temporary solution to make score_x_i to overweight other scores to get significant clusters
        score_x_i = -np.log(score_x_i)
        # CHECK in evalue section how to initialize this score
        strand = str(mapped_results.iloc[i,1])
        gap = abs(int(mapped_results.iloc[i,2]) - int(mapped_results.iloc[i-1,2])) - 1
        print('gap =', gap)
        # CHECK if that's strand of the previous gene or next
        # FIX with OR statement (re.findall('\+|\-', test))
        if ''.join(re.findall('\+', strand))==init_strand:
            f_strand_flip = 0
            print('same strand')
        else:
            if ''.join(re.findall('\-', strand))==init_strand:
                f_strand_flip = 0
                print('same strand')
            else:
                f_strand_flip = 1
        
        # updating previous gene strand (current gene = previous for next for loop iter)
        if re.findall('\+', strand)==[]:
            init_strand= '-'
        if re.findall('\-', strand)==[]:
            init_strand= '+'
        
        print('scores')
        print(score_i_minus_1_cluster)
        print(f_strand_flip*d_strand_flip_penalty)
        print(score_x_i)
        print('  ')
        print(score_i_minus_1_cluster - f_strand_flip*d_strand_flip_penalty + score_x_i)
        print(max(0, score_x_i))
        # CHECK gap penalty
        if (score_i_minus_1_cluster - f_strand_flip*d_strand_flip_penalty - gap_penalty*gap + score_x_i) > max(0, score_x_i):
            score_i_cluster = score_i_minus_1_cluster - f_strand_flip*d_strand_flip_penalty - gap_penalty*gap + score_x_i
            print('proceed', score_i_cluster, score_max_cluster)
            print('score i-1', score_i_minus_1_cluster)
            
            if score_i_cluster > score_max_cluster:
                print('first')
                score_max_cluster = score_i_cluster
                # CHECK if correct, CHANGE to look better
                #i_1_cluster_end = str(mapped_results.iloc[i,1])
                i_1_cluster_end = str(mapped_results.iloc[i,1])[str(mapped_results.iloc[i,1]).find(''.join(re.findall('\+|\-', str(mapped_results.iloc[i,1]))))+1:]
                # CHECK changing of score_i_minus_1_cluster
                score_i_minus_1_cluster = score_i_cluster
        else:
            print('second')
            score_i_cluster = score_x_i
            # CHECK if correct, CHANGE to get right coord
            # CHECK, maybe I have to take i-1 coord? otherwise
            # it seems like I take the left coord of 1st non-cluster
            # prot - ASK Johannes
            
            if score_max_cluster > score_min_cluster:
                print('1st append')
                cluster_matches.append((i_0_cluster_start,
                i_1_cluster_end, score_max_cluster))
                score_max_cluster = 0
            # CHECK if ok to shift it here from above if 
            # CHANGE to more well-looking finding of a strand/coord
            i_0_cluster_start = str(mapped_results.iloc[i,1])[:str(mapped_results.iloc[i,1]).find(''.join(re.findall('\+|\-', str(mapped_results.iloc[i,1]))))]
        print('max and min scores', score_max_cluster, score_min_cluster)
        print('cluster coord', i_0_cluster_start, i_1_cluster_end)
        print(cluster_matches)

    if score_max_cluster > score_min_cluster:
        print('2nd append')
        cluster_matches.append((i_0_cluster_start,
         i_1_cluster_end, score_max_cluster))
    # add more to cluster matches table? prot id?
    print(cluster_matches)
    # ADD return of the changed ResultsMapping object? (with added scores?)
    return cluster_matches


# it re-defines the scores which I have got in the first iteration 
# CHECK if i have to read some combined\changed file of ResultsMapping
# ADD s0 
def update_scores_for_cluster_matches(cluster_matches):
    significant_clusters = cluster_matches 
    # ADD query id to mapped results
    mapped_results = mapped_res.res_map_to_header
    results = mapped_res.search_result_file
    # iterate through query ids in mapped results (ADD this column)
    # CHECK if these are right columns
    L = len(significant_clusters)
    print(L)
    bias = 0
    for query_id in mapped_results.iloc[:,3]:
        print(query_id)
        m_x = mapped_results.iloc[:,3][mapped_results.iloc[:,3] == query_id].shape[0]
        M_x = results.iloc[:,0][results.iloc[:,0] == query_id].shape[0]
        # CHECK if true
        l = m_x
        print(m_x, M_x, l)
        score_x = np.log(np.divide(np.divide(m_x, l), np.divide(M_x, L))) - bias
    pass


def set_strand_flip_penalty():

    pass


def update_query_profiles():
    pass


def add_new_proteins():
    pass


def main():
    #make_profiles()

    # CHECK if it works with multihitdb (just from command line it worked)
    # CHECK why there are more results with multihitdb (target is converted to profiles??)
    #run_search()

    # this class is to have order and strand for target proteins
    # FIX best query hit is needed
    #mapped_res = ResultsMapping.map_target_to_coord()
    #print(mapped_res.res_map_to_header)

    #cluster_matches = find_clusters()
    # CHECK if it is optimal to divide scores update to few functions
    #update_scores_for_cluster_matches(cluster_matches)
    #set_strand_flip_penalty()
    #update_query_profiles()
    #add_new_proteins()
    pass


if __name__ == "__main__":

    print("starting")

    print('type number of iterations')
    #iterations = int(input())
    iterations = 1

    files = FilePath.get_path()
    print(files.query_db)
    print(files.target_db)
    print(files.res)

    # FIX to be the right order of functions (should be after run_search())
    mapped_res = ResultsMapping.map_target_to_coord()

    while iterations > 0:
        main()
        iterations -= 1


