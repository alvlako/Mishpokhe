import numpy as np
from io import StringIO
import pandas as pd

import os
import subprocess

# version 1 to get enrichment (whether proteins are enriched in clusters, not alone)
# and add these protein to the query for the searches
# CHANGE to take unordered dataset

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


class SearchRes():
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
    

# Currently these steps are run for target against target
def run_clustersearch():
    print('running multihitsearch')
    subprocess.call(['mmseqs', 'multihitsearch', cls_files.target_db,
    cls_files.target_db, str(cls_files.res + str("_cls")), 'tmp'])


# Currently these steps are run for target against target
# ADD different search modes, currentle - prot against prot
def run_search():
    print('running mmseqs search')
    subprocess.call(['mmseqs', 'search', norm_files.target_db,
    norm_files.target_db, str(norm_files.res + str("_norm")), 'tmp', '--search-type', '3'])


def count_matches_proportion():
    # transfer in main? because 1 is only for 1st iter (what?)
    enrichment_score = 1
    os.chmod('/Users/Sasha/Desktop/bash_mishpokhe.sh', 0o755)
    cls_res_file = str(cls_files.res + str("_cls"))
    norm_res_file = str(norm_files.res + str("_norm"))
    subprocess.call(["zsh", "/Users/Sasha/Desktop/bash_mishpokhe.sh",
     cls_res_file, norm_res_file, cls_files.res])


# Maybe I do not get how the scores are calculated (for bash script as well)
# I feel it when read the formula 2 for calculating the score
# for proteins with no match
# for the time just ignore 0 and no matched proteins FIX
# THINK is it right that clustersearch results and normal search results have the same ids?
# and i can just merge the results by id?

# calculate score from clustersearch \ normal search matches
# according to the formula
def calculate_score():
    scores_file = str(cls_files.res + str("_all_score"))
    scores = np.loadtxt(scores_file, delimiter="\t")
    print(scores)
    print("Type b bias")
    # Is it okay that it is a float? Make checks for correct unput?

    #bias = float(input())
    bias = float(0.00001)
    # ignoring zeros for now FIX
    with np.errstate(divide='ignore'):
        enrich_score = np.log(np.divide(scores[:,2],scores[:,3])) - bias
    enrich_scores_arr = np.append(scores, enrich_score[:,None], axis=1)
    print(enrich_scores_arr)
    # get only lines with the score above something
    # to extract cluster-enriched proteins
    # make this threshold learning from the data?
    #score_threshold = np.mean(enrich_scores_arr[:,4]) # too many prot remain. why?
    score_threshold = 0.09
    bool_mask = (enrich_scores_arr[:, 4] > score_threshold)
    tmp_enrich = enrich_scores_arr[bool_mask, :]
    enriched_prot_ids = tmp_enrich[:,0].astype(int)
    print(enriched_prot_ids)
    print(enrich_scores_arr[:,1].size)
    print(enriched_prot_ids.size)
    np.savetxt(str(cls_files.res + str("_enrich_prot_1iter")), enriched_prot_ids,
     delimiter='\t', fmt='%i')


def update_query():
    subprocess.call(['mmseqs', 'search', norm_files.query_db,
    norm_files.target_db, str(norm_files.res + str("_norm")), 'tmp', '--search-type', '3'])
    pass


# Need to be changed, saving and loading files again is not optimal at all
# At this try I will take the prots from multihitdb
# I thought I had to map the protein ids from results to lookup
# and find an index of the line in lookup and extract the corresponding line
# in the db, but I dont need the lookup, the results id is already an index 
# of the sequence line
# WHY db.shape = 165 while lookup is 164? Last line in db is empty?
# subset test multihit db to get only enriched proteins - do subset for
# the db file, lookup, index and _h
# What is wrong with the last line in _h? Consists 1 col instead of 2
# Rewrite accordingly to the new scheme
# Is this step for query against target approach? THINK
def add_new_proteins():
    new_prot_id = np.loadtxt(str(cls_files.res + str("_enrich_prot_1iter")))
    target_db = np.genfromtxt(str(cls_files.target_db), dtype = None, delimiter="\t")
    target_db_lookup = np.genfromtxt(str(cls_files.target_db)+str(".lookup"), dtype = None, delimiter="\t", encoding=None)
    target_db_index_tmp = pd.read_csv(str(cls_files.target_db)+str(".index"), sep='\t',header=None)
    # Is it ok to ignore raised errors to ignore problem with the last line??
    target_db_h_tmp = pd.read_csv(str(cls_files.target_db)+str("_h"), sep='\t',header=None)
    
    target_db_index = pd.melt(target_db_index_tmp).to_records()
    target_db_h = pd.melt(target_db_h_tmp).to_records()
    
    # SOMETHING is weird with the _h FILE FIX!!!
    target_db_enrich_prot_db = np.take(target_db, new_prot_id.astype(int))
    target_db_enrich_prot_lookup = np.take(target_db_lookup, new_prot_id.astype(int))
    target_db_enrich_prot_index = np.take(target_db_index, new_prot_id.astype(int))
    target_db_enrich_prot_h = np.take(target_db_h, new_prot_id.astype(int))

    np.savetxt(str(cls_files.res + str("_enrich_prot_db")), target_db_enrich_prot_db.astype(np.str_), fmt='%s')
    np.savetxt(str(cls_files.res + str("_enrich_prot_db.lookup")), target_db_enrich_prot_lookup, delimiter="\t", fmt='%s')
    np.savetxt(str(cls_files.res + str("_enrich_prot_db.index")), target_db_enrich_prot_index, delimiter="\t", fmt='%s')
    np.savetxt(str(cls_files.res + str("_enrich_prot_db_h")), target_db_enrich_prot_h, delimiter="\t", fmt='%s')
    #make for multihitdb. Do i need this merge if that's target against target?
    # CHECK might not work for other steps\searches because the results db is actually
    # a part of the real db, not the normal db with all indices
    subprocess.call(['mmseqs', 'mergedbs', cls_files.target_db,
     str(cls_files.target_db+ str("_added")), str(cls_files.res + str("_enrich_prot_db"))])


# That's for target against target approach
# improve SPEED
def get_enrich_cluster_target_prot():
    f1 = open(str(cls_files.res + str("_enrich_prot_db"))) 
    f2 = open(str(cls_files.res + str("_enrich_prot_db.lookup")))
    f3 = open(str(cls_files.res) + str("_enrich_prot.faa"), 'w')
    l1 = f1.readlines()
    l2 = f2.readlines()
    for i1,i2 in zip(l2, l1):
        f3.writelines(str('>') + i1)
        f3.writelines(i2)
    f1.close()
    f2.close()
    f3.close()
    #target_db = np.genfromtxt(str(cls_files.target_db), dtype = None, delimiter="\t")
    #target_db_lookup = np.genfromtxt(str(cls_files.target_db)+str(".lookup"), dtype = None, delimiter="\t", encoding=None)
    #pass


# Should i remove old results each step? THINK
# FIX!!! No k-mer could be extracted for the database py_res_enrich_prot.faa_db.Maybe the sequences length is less than 14 residues
def unordered_against_enriched():
    #subprocess.call(['mmseqs', 'createdb', str(cls_files.res + str("_enrich_prot.faa")),
     #str(cls_files.res + str("_enrich_prot.faa_db"))])
    #subprocess.call(['mmseqs', 'search', norm_files.query_db,
     #str(cls_files.res + str("_enrich_prot.faa_db")), str(norm_files.res + str("_ag_enriched_target_res")), 'tmp'])
    subprocess.call(['mmseqs', 'createdb', 'entest.faa',
     'entest.faa_db'])
    subprocess.call(['mmseqs', 'search', norm_files.query_db,
     'entest.faa_db', 'query_ag_entest_res', 'tmp'])
    #pass


def main():
    # Part 1 to find target ordered genomes proteins enriched in clusters
    # add proteins REDO to add to the query unordered set, not to the target
    # Currently it is target against target

    #run_clustersearch()
    #run_search()
    #count_matches_proportion()

    # SCORE calculation did not work! FIX!!!

    #calculate_score()

    #get_enrich_cluster_target_prot()

    # Is this step for query against target approach? Well, still useful to get proteins as 
    #it constructs the db of the enriched proteins

    #add_new_proteins()

    # Part 2 to find query proteins matches to target cluster-enriched proteins
    # and to add the matches to the query

    unordered_against_enriched()
    #update_query()

    #pass


# CHANGE to ask only for 1 result file
# CHANGE for normal iterations scheme
if __name__ == "__main__":

    print("starting")

    print('type number of iterations')
    #iterations = int(input())
    iterations = 1

    # it asks twice for the file paths, for multihitdb for 
    # the clustersearch and normal db for the normal search
    # that's for TARGET all against all
    # query is the same, CHANGE to not ask about it
    cls_files = FilePath.get_path()
    print(cls_files)
    print(cls_files.query_db)
    print(cls_files.target_db)
    print(cls_files.res)
    norm_files = FilePath.get_path()

    while iterations > 0:
        main()
        iterations -= 1
    
    # it asks for QUERY unordered proteins set db. Not needed if that's query
    # and the initial steps are target against target?
    
    #unord_files = FilePath.get_path()
    #print(unord_files)
    #print(unord_files.query_db)
    #print(unord_files.target_db)
    #print(unord_files.res)


