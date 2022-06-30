import numpy as np
import pandas as pd

import os
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
def run_search():
    print('running mmseqs profile search')
    subprocess.call(['mmseqs', 'search', files.target_db,
    files.query_db + '_clu' + '_rep' + '_profile', files.res + '_prof_search',
     'tmp', '-a'])


def find_clusters():
    pass


def update_scores():
    pass


def update_query_profiles():
    pass


def add_new_proteins():
    pass


def main():
    make_profiles()
    run_search()
    #find_clusters()
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

    while iterations > 0:
        main()
        iterations -= 1


