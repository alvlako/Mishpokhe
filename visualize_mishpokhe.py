from cmath import inf
import ctypes as ct
from ctypes import *
import numba as nb
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.figure_factory as ff
import plotly.io as pio

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
import packaging
import networkx as nx
import matplotlib.pyplot as plt

path_to = 'res_terpene_acidobacteria2iter_res_3_iter_sign_clusters_enrich_stat_filtered_clu_filter'
path_to_test = path_to
# ast.literal_eval is used to read rows looking like [1,2,3] as python list
clusters_stat = pd.read_csv(path_to_test, dtype=None, sep='\t')
print(clusters_stat)
clusters_stat = pd.read_csv(path_to_test, dtype=None, sep='\t',
       converters={'query_prots':ast.literal_eval,
      'coordS1':ast.literal_eval, 'coordS2':ast.literal_eval, 'target_prots':ast.literal_eval})
# '
# 'strand':ast.literal_eval, 


clusters_stat['targer_string'] = [','.join(map(str, l)) for l in clusters_stat['target_prots']]
clusters_stat['target'], clusters_stat['rest'] = clusters_stat['targer_string'].str.split(',', n=1).str
clusters_stat['target_genome'], clusters_stat['rest2'] = clusters_stat['target'].str.rsplit('_', n=1).str

short = clusters_stat[["target_genome", "target_prots", "coordS1", "coordS2", "query_prots"]]
short["cluster_id"] = short["target_genome"]

# to have cluster ids
curr_genome = ''
prev_genome = ''
n = 1
for i in range(0, len(short.iloc[:, 1])):
    curr_genome = short.loc[i, 'target_genome']
    #print('curr_genome', curr_genome, 'prev_genome', prev_genome)
    if curr_genome == prev_genome:
       short.loc[i, 'cluster_id'] = short.loc[i, 'target_genome'] + f'_clu_{n}'
    else:
       n = 1
       short.loc[i, 'cluster_id'] = short.loc[i, 'target_genome'] + f'_clu_{n}'
    n = n + 1
    prev_genome = curr_genome

pd.set_option('display.max_columns', None)
print(short)

short_exp = short.explode(['target_prots', 'coordS1', 'coordS2', 'query_prots']).reset_index(drop=True)

# to have shifted protein coords per cluster, not per genome
clu_start = 0
short_exp["coord1_shift"] = 0
short_exp["coord2_shift"] = 0
short_exp["coord_diff"] =  short_exp["coordS2"]  - short_exp["coordS1"] 
curr_clu = ''
prev_clu = ''
for i in range(0, len(short_exp.iloc[:, 1])):
    curr_clu = short_exp.loc[i, 'cluster_id']
    #print('curr_genome', curr_genome, 'prev_genome', prev_genome)
    if curr_clu == prev_clu:
       short_exp.loc[i, "coord1_shift"] = short_exp.loc[i, "coordS1"] - clu_start
       short_exp.loc[i, "coord2_shift"] = short_exp.loc[i, "coordS2"] - clu_start
    else:
       short_exp.loc[i, "coord2_shift"] = short_exp.loc[i, "coord_diff"] 
       clu_start = short_exp.loc[i, "coordS1"]
    prev_clu = curr_clu

print(short_exp)


short_exp["target_genome_intern_id"] = short_exp.groupby('target_genome', sort=False).ngroup() + 1


print(short_exp)

pio.renderers.default = "chrome"

df = [dict(Task="genome1", Start=20, Finish=100, Resource='Complete'),
      dict(Task="genome1", Start=120, Finish=400, Resource='Incomplete'),
      dict(Task="genome2", Start=20, Finish=80, Resource='Not Started'),
      dict(Task="genome2", Start=140, Finish=200, Resource='Complete'),
      dict(Task="genome3", Start=20, Finish=200, Resource='Not Started'),
      dict(Task="genome3", Start=200, Finish=250, Resource='Not Started'),
      dict(Task="genome3", Start=300, Finish=420, Resource='Not Started'),
      dict(Task="genome4", Start=20, Finish=150, Resource='Complete')]

#colors = {'Not Started': 'rgb(220, 0, 0)','Incomplete': (1, 0.9, 0.16),'Complete': 'rgb(0, 255, 100)'}

#fig = ff.create_gantt(df, colors=colors, index_col='Resource', show_colorbar=True,group_tasks=True)
#fig.show()

clusters_list = []
for i in range(0, len(short_exp.iloc[:, 1])):
    genome = short_exp.loc[i, 'cluster_id']
    prot = short_exp.loc[i, 'target_prots']
    query_prot = short_exp.loc[i, 'query_prots']
    coord1 = short_exp.loc[i, 'coord1_shift']
    coord2 = short_exp.loc[i, 'coord2_shift']
    clusters_list.append({'Task': genome, 'Start': coord1, 'Finish': coord2, 'Resource': query_prot})

#colors = {'Not Started': 'rgb(220, 0, 0)', 'Incomplete': (1, 0.9, 0.16), 'Complete': 'rgb(0, 255, 100)'}
#colors=colors,

colors = dict()
c = 0
for i in short_exp["query_prots"].unique():
    colors[i] = px.colors.qualitative.Dark24[c]
    c = c + 1
    if c > 23:
        c = 0


#print(clusters_list)

def create_gantt_clusters():
    fig = ff.create_gantt(clusters_list, colors = colors, index_col='Resource',   show_colorbar=True, group_tasks=True, title='Mishpokhe clusters')
    fig.update_layout(xaxis_type='linear')
    #fig.layout.title = None
    #fig.layout.xaxis.rangeselector = None
    fig.show()

create_gantt_clusters()
print(x)


import plotly.graph_objects as go
fig = go.Figure()

genomes_number = len(short_exp["target_genome_intern_id"].unique())

fig.add_trace(go.Scatter(x=[short_exp["coordS2"].max()+300 for i in range(0, genomes_number)], y=[i for i in range(0,genomes_number)], mode="lines", opacity=0))
for i in range(0, genomes_number):
    fig.add_hline(y=i)

#fig.show()
#print(x)

opt_rect_height = genomes_number / (short_exp["coordS2"].max()/short_exp["coord_diff"].min())
opt_rect_height = 0.5
for x0, y0, x1, y1 in zip(short_exp["coordS1"].tolist(), [i-opt_rect_height for i in short_exp["target_genome_intern_id"].tolist()], short_exp["coordS2"].tolist(), [i+opt_rect_height for i in short_exp["target_genome_intern_id"].tolist()]):
    #print(x0, y0, x1, y1 )
    fig.add_shape(type="rect",
        x0=x0, y0=y0, x1=x1, y1=y1,
        fillcolor="green",
        opacity=0.7,
        layer="below",
        line_width=0
    )
fig.update_shapes(dict(xref='x', yref='y'))
fig.show()

