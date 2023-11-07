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

short = clusters_stat[["target_genome", "target_prots", "coordS1", "coordS2", "query_prots", "arc_clu_id"]]
short['list_new_scoreS_enrich'] = clusters_stat['list_new_scoreS_enrich'].apply(lambda x: x.split(',')[1:])
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
print('SHORT', short)

# for some reason there are spare elements (probably 0 at the beginning, have to figure out)
# I have to fix the initial code
short.loc[i, 'list_new_scoreS_enrich_fixed'] = ""
for i in range(0, len(short.iloc[:, 1])):
    if len(short.loc[i, 'coordS1']) != len(short.loc[i, 'list_new_scoreS_enrich']):
        if len(short.loc[i, 'coordS1']) < len(short.loc[i, 'list_new_scoreS_enrich']):
            short.at[i, 'list_new_scoreS_enrich_fixed']  = short.loc[i, 'list_new_scoreS_enrich'][1:]
        else:
            short.at[i, 'list_new_scoreS_enrich_fixed']  = short.loc[i, 'list_new_scoreS_enrich']
            for k in range(0, len(short.loc[i, 'coordS1']) - len(short.loc[i, 'list_new_scoreS_enrich'])):
                short.at[i, 'list_new_scoreS_enrich_fixed'].append(0)
    else:
        short.at[i, 'list_new_scoreS_enrich_fixed']  = short.loc[i, 'list_new_scoreS_enrich']

for i in range(0, len(short.iloc[:, 1])):
    if len(short.loc[i, 'coordS1']) != len(short.loc[i, 'list_new_scoreS_enrich_fixed']):
        print(short.loc[i, 'coordS1'])
        print(short.loc[i, 'list_new_scoreS_enrich_fixed'])
        short.at[i, 'list_new_scoreS_enrich_fixed'] = short.loc[i, 'list_new_scoreS_enrich_fixed'][1:]

short_exp = short.explode(['target_prots', 'coordS1', 'coordS2', 'query_prots','list_new_scoreS_enrich_fixed']).reset_index(drop=True)


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
short_exp["cluster_id_intern_id"] = short_exp.groupby('cluster_id', sort=False).ngroup() + 1


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
annots_list = []
clu_number = len(short_exp['cluster_id_intern_id'].unique())+1 
#short_exp = short_exp.dropna()
short_exp = short_exp.sort_values(['arc_clu_id', 'cluster_id']).reset_index(drop=True)

#pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
print(short_exp)
#print(x)

prev_genome = ''
curr_clu_index = 0
for i in range(0, len(short_exp.iloc[:, 1])):
    genome = short_exp.loc[i, 'cluster_id']
    prot = short_exp.loc[i, 'target_prots']
    query_prot = short_exp.loc[i, 'query_prots']
    coord1 = short_exp.loc[i, 'coord1_shift']
    coord2 = short_exp.loc[i, 'coord2_shift']
    arc_clu_id = short_exp.loc[i, 'arc_clu_id']
    clusters_list.append({'Task': genome, 'Start': coord1, 'Finish': coord2, 'Resource': query_prot, 'Description': prot, "Ind": arc_clu_id})

    x1 = coord1 + ((coord2 - coord1) / 2)
    print('coord1', coord1, 'coord2', coord2, 'x1', x1)
    #curr_clu_index
    if genome != prev_genome:
        print('genome', genome, 'prev_genome', prev_genome)
        prev_genome = genome
        curr_clu_index = curr_clu_index + 1
    # that is to reverse the order
    y1 = clu_number - curr_clu_index - 1
    #y1 = clu_number - n
    try:
        text1 = round(float(short_exp.loc[i, 'list_new_scoreS_enrich_fixed']), 1)
        #text1 = f'curr_clu_index is {curr_clu_index}'
    except ValueError:
        text1 = ''
    annots_list.append({'x': x1, 'y': y1, 'text': text1, 'showarrow': False, 'font': {'color': 'black'}})
    

#colors = {'Not Started': 'rgb(220, 0, 0)', 'Incomplete': (1, 0.9, 0.16), 'Complete': 'rgb(0, 255, 100)'}
#colors=colors,

colors = dict()
c = 0
short_exp["query_colors"] = ''
for i in short_exp["query_prots"].unique():
    #colors[i] = f'rgba{px.colors.hex_to_rgb(px.colors.qualitative.Dark24[c]) + (1,)}'
    colors[i] = px.colors.qualitative.Dark24[c]
    print(colors[i])
    c = c + 1
    if c > 23:
        c = 0

for i in range(0, len(short_exp.iloc[:, 1])):
    short_exp.loc[i, "query_colors"] = colors[short_exp.loc[i, "query_prots"]]

#short_exp['opacity'] = (1 - (1/pd.to_numeric(short_exp['list_new_scoreS_enrich_fixed']))).replace(-np.Inf, 0)
#short_exp['opacity'][short_exp['opacity'] < 0] = 0
#print(short_exp['opacity'])

#print(clusters_list)

def create_gantt_clusters():

    fig = ff.create_gantt(clusters_list, colors = colors, index_col='Resource',   show_colorbar=True, group_tasks=True, title='Mishpokhe clusters', showgrid_x=True, showgrid_y=True)
    fig.update_layout(xaxis_type='linear')
    #fig.update_traces(opacity=0.1, selector={'fill':'toself'})
    #fig.layout.title = None
    #fig.layout.xaxis.rangeselector = None
    my_fig = fig.to_dict()
    

    fig['layout']['annotations'] = annots_list
    config = {'scrollZoom': True}
    #fig.update_layout(xaxis=dict(rangeslider=dict(visible=True),type="linear"))
    #pio.show(my_fig)
    fig.show(config=config)

create_gantt_clusters()
print(x)


import plotly.graph_objects as go
fig = go.Figure()

clusters_number = len(short_exp["cluster_id"].unique())

fig.add_trace(go.Scatter(x=[short_exp["coord2_shift"].max()+300 for i in range(0, clusters_number)], y=[i for i in range(0,clusters_number)], mode="lines", opacity=0))
for i in range(0, clusters_number):
    fig.add_hline(y=i)

#fig.show()
#print(x)

opt_rect_height = clusters_number / (short_exp["coordS2"].max()/short_exp["coord_diff"].min())
opt_rect_height = 0.2
for x0, y0, x1, y1, op, col, nam in zip(short_exp["coord1_shift"].tolist(), [i-opt_rect_height for i in short_exp["cluster_id_intern_id"]], short_exp["coord2_shift"].tolist(), [i+opt_rect_height for i in short_exp["cluster_id_intern_id"]], short_exp['opacity'], short_exp['query_colors'], short_exp['query_prots'] ):
    #print(x0, y0, x1, y1 )
    fig.add_shape(type="rect",
        x0=x0, y0=y0, x1=x1, y1=y1,
        fillcolor=col,
        opacity=op,
        layer="below",
        line_width=0,
        name = nam,
        showlegend=True
    )
fig.update_shapes(dict(xref='x', yref='y'))
fig.show()

