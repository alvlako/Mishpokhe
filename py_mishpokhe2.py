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
    # FIX something wrong local id (4294967295) >= db size (4)
    #subprocess.call(['mmseqs', 'convertalis', files.query_db + '_clu' + '_rep' + '_profile',
     #files.target_db, files.res + '_prof_search',
      #files.res + '_prof_search' +'.m8'])


# The n of prots in db matches the n of gff CDS CDS!!! entries
# FIX best query hit is needed
# CHECK if the mapping correct 
# CHECK if it quaranteed that indices in results file are same order as in lookup (also in order) 
# AND in _h and the order is always the same
# FIX !!!!!! some real ids are read as NaN
class ResultsMapping:
     """
     Mapping results to their coordinates and strands

     """
     def __init__(self, search_result_file, target_db_lookup, target_db_h, res_map_to_header, ind_list):
        self.search_result_file = search_result_file
        self.target_db_lookup = target_db_lookup
        self.target_db_h = target_db_h
        self.res_map_to_header = res_map_to_header
        self.ind_list = ind_list


     @classmethod
     def map_target_to_coord(self):
        #search_result_file = pd.read_csv(str(files.res + '_prof_search'), dtype={'int':'float'}, sep='\t')
        search_result_file = pd.read_csv('/Users/Sasha/Documents/GitHub/mishpokhe_test/py2_multihit_res', dtype={'int':'float'}, sep='\t')
        # given the target db is multihitdb. Do i need lookup?
        target_db_lookup = np.genfromtxt(str(files.target_db)+str(".lookup"), dtype = None, delimiter="\t", encoding=None)
        target_db_h = pd.read_csv(str(files.target_db)+str("_h"), sep='\t',header=None)
        #print(search_result_file.iloc[:, 0])
        #print(target_db_h)
        ind_list = search_result_file.iloc[:, 0].dropna().astype(int)
        #print(ind_list)
        res_map_to_header = target_db_h.iloc[ind_list]
        #print(type(ind_list))
        res_map_to_header['ind'] = ind_list.values
        #print(res_map_to_header)
        return self(search_result_file, target_db_lookup, target_db_h, res_map_to_header, ind_list)



# ?? is it correct that L is set not of all the target proteins
# but only those which were matched and got to results???
# FIX !!!!!! some real ids are read as NaN
# FIX ! to process not the initial results but the best hit results
# CHECK whether you need to change the scores
def find_clusters():

    mapped_results = mapped_res.res_map_to_header
    results = mapped_res.search_result_file
    index_list = mapped_res.ind_list
    print(results)
    print(mapped_results)
    print(index_list)

    #Algorithm 1 - iterate through target prot
    print(results.iloc[:, 0].size)

    cluster_matches = list()
    # CHECK if score max cluster set up correct
    score_max_cluster = 0
    # OPTIMIZE how to set the strand penalty
    d_strand_flip_penalty = 0.5
    score_min_cluster = 0
    score_i_minus_1_cluster = 0
    # CHECK if correct
    strand = str(mapped_results.iloc[0, 2])
    if re.findall('\+', strand)==[]:
        init_strand= '-'
    if re.findall('\-', strand)==[]:
        init_strand= '+'
    print(init_strand)

    print(len(results.iloc[:, 0]))
    print(results.iloc[:, 0])
    for i in range(0, len(results.iloc[:, 0])-1):
        print(i)
        print(results.iloc[i, :])
        score_x_i = float(results.iloc[i, 3])
        # CHECK in evalue section how to initialize this score
        strand = str(mapped_results.iloc[i, 2])
        # CHECK if that's strand of the previous gene or next
        if re.findall('\+', strand)==init_strand:
            f_strand_flip = 0
        if re.findall('\-', strand)==init_strand:
            f_strand_flip = 0
        else:
            f_strand_flip = 1
        
        if (score_i_minus_1_cluster - f_strand_flip*d_strand_flip_penalty + score_x_i) > max(0, score_x_i):
            score_i_cluster = score_i_minus_1_cluster - f_strand_flip*d_strand_flip_penalty + score_x_i
            if score_i_cluster > score_max_cluster:
                score_i_cluster = score_i_cluster
                i_1_cluster_end = i
            else:
                score_i_cluster = score_x_i
                i_0_cluster_start = i
                if score_max_cluster > score_min_cluster:
                    cluster_matches.append((i_0_cluster_start,
                    i_1_cluster_end, score_max_cluster))
                    score_max_cluster = 0
    if score_max_cluster > score_min_cluster:
        cluster_matches.append((i_0_cluster_start,
         i_1_cluster_end, score_max_cluster))
    print(cluster_matches)
    return cluster_matches


def update_scores():
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

    cluster_matches = find_clusters()
    #update_scores()
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


