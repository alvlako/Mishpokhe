# That is a script that should generate figures for mishpokhe paper

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
import random
import re
import subprocess
import sys
from sys import stdout
import time
import packaging
import networkx as nx
import matplotlib.pyplot as plt
import squarify
import seaborn as sns


# Figure 1
# That is a figure that should visualize via "clustering" (not a real, just for visualization) how many defence systems out of all (including new not only the seed ones) i have found

# Let's visualize first all the defence systems found in the data by padloc
# I have to generate random coordinates for them

# Let's load the padloc ground truth 
path_to_padloc = '/Users/Sasha/Downloads/padloc_res_ncbi_bact_chr_repr_short'
padloc_systems = pd.read_csv(path_to_padloc, dtype=None, sep='\t', header = None)
padloc_systems.columns = ['genome', 'coord1', 'coord2', 'type']
print(padloc_systems)
uniq_padloc_systems = pd.DataFrame(padloc_systems['type'].value_counts()).reset_index()
uniq_padloc_systems.columns = ['type', 'counts']
# make a column for colors
uniq_padloc_systems['colors'] = uniq_padloc_systems.index
# Let's assign intervals for random generation
number_of_systems = len(uniq_padloc_systems)
max_number_per_system = int(uniq_padloc_systems['counts'].max())
print('max_number_per_system', max_number_per_system)
multiplier = max_number_per_system
width = number_of_systems*multiplier
height = number_of_systems*multiplier
step_x = int(width/number_of_systems)
step_y = int(height/number_of_systems)
print('width, height, step_x, step_y', width, height, step_x, step_y)
print('number_of_systems', number_of_systems)
uniq_padloc_systems['x_interval1'] = [i for i in range(0, width,step_x)]
uniq_padloc_systems['x_interval2'] = uniq_padloc_systems['x_interval1'] + step_x
list_y_1 = [i for i in range(0, height,step_y)]
random.shuffle(list_y_1)
uniq_padloc_systems['y_interval1'] = list_y_1
uniq_padloc_systems['y_interval2'] = uniq_padloc_systems['y_interval1'] + step_y
uniq_padloc_systems['x_coord_array'] = ''
uniq_padloc_systems['y_coord_array'] = ''
for i in range(0,number_of_systems):
    startx = int(uniq_padloc_systems.iloc[i]['x_interval1'])
    stopx = int(uniq_padloc_systems.iloc[i]['x_interval2'])
    starty = int(uniq_padloc_systems.iloc[i]['y_interval1'])
    stopy = int(uniq_padloc_systems.iloc[i]['y_interval2'])
    number_of_randoms = int((uniq_padloc_systems.iloc[i]['counts']))
    print('startx, stopx, starty, stopy, number_of_randoms', startx, stopx,starty, stopy, number_of_randoms)
    list_of_randomsx = list(np.random.randint(low=startx, high=stopx, size=(number_of_randoms,)))
    list_of_randomsy = list(np.random.randint(low=starty, high=stopy, size=(number_of_randoms,)))
    uniq_padloc_systems.iat[i, uniq_padloc_systems.columns.get_loc('x_coord_array')] = list_of_randomsx
    uniq_padloc_systems.iat[i, uniq_padloc_systems.columns.get_loc('y_coord_array')] = list_of_randomsy

print(uniq_padloc_systems)
# Now i have to explode the column with coordinates so I have coordinates for every occurence of the system
padloc_systems_exploded = uniq_padloc_systems.explode(['x_coord_array','y_coord_array']).reset_index()
print('padloc_systems_exploded', padloc_systems_exploded)

#plt.scatter(data = padloc_systems_exploded, x='x_coord_array', y='y_coord_array', c='colors', alpha=0.5, label='type')
#plt.legend()
#plt.show()

# plot
# let's try to first make a column that has a more broad classification of types
padloc_systems_exploded['type_to_split'] = padloc_systems_exploded['type'].astype(str)+'_'
print(padloc_systems_exploded)
print(padloc_systems_exploded['type_to_split'].str.split('_',expand=True)[0])
padloc_systems_exploded['general_type'] = padloc_systems_exploded['type_to_split'].str.split('_',expand=True)[0]

pio.renderers.default = "chrome"

fig = px.scatter(padloc_systems_exploded, x="x_coord_array", y="y_coord_array", color="general_type", size='counts')
fig.show()

