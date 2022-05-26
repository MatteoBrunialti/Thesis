import csv
import os
import statistics as sta
from Bio import Phylo
import io
import itertools
import numpy as np
import pandas as pd
from itertools import product
from itertools import permutations

maindir = '/data/matteo/OrthoFinder/Results_May20'
file_tsv = f"{maindir}/Orthogroups/Orthogroups.tsv"
count_dir = f'{maindir}/Orthogroups/Orthogroups.GeneCount.tsv'

def read(file_tsv, count_dir):
    df = pd.read_table(file_tsv)
    finaldict = {}
    column_no_dot = [x.replace('.', '_') for x in df.columns]
    num = pd.read_table(count)
    file_num = num[num['Total'] == 3].index[0]

    #creating the base dictionary
    for i in [f'{column_no_dot[x]}|{column_no_dot[y]}' for x in range(1, len(column_no_dot)) for y in range(x + 1, len(column_no_dot))]:
        finaldict[i] = []


    for i in range(file_num):  # TO DO THE RANGE CHECK FOR FILE
        OG = df.iloc[i,0]
        small_store = []
        for j in range(1, len(df.columns)):
            if type(df.iloc[i,j]) == str:
                small_store.append([f'{column_no_dot[j]}_{x}'.replace(' ', '') for x in df.iloc[i, j].split(',')])

        tree = Phylo.read(f'{maindir}/Resolved_Gene_Trees/{OG}_tree.txt', "newick")
        if len(small_store) > 1:
            for one in range(len(small_store)):
                for two in range(one + 1, len(small_store)):
                    # here modify so to select species vs species + distance
                    var = [(x, y, tree.distance(x, y)) for x in small_store[one] for y in small_store[two]]
                    # joining the two names of the tuple to match the dictionary
                    obj3 = '_'.join(var[0][0].split('_')[:-1]) + '|' + '_'.join(var[0][1].split('_')[:-1])
                    finaldict[obj3] += var
    return finaldict

def statistics(dictionary, verbosity = 0):
    outliers = {}
    for j in dictionary:  # j = the key; dictionary[k]  = the value
        mavalue = [dictionary[j][x][2] for x in range(len(dictionary[j]))]
        left_arm = np.percentile(mavalue, 1)
        outliers[j] = [dictionary[j][x] for x in range(len(mavalue)) if mavalue[x] <= left_arm]
        if verbosity == 0:
            return outliers
        elif verbosity == 1:
            print(f'species: {j}')
            print(f'mean: {sta.mean(mavalue)}')
            print(f'mode:{sta.mode(mavalue)}')
            print(f'median: {sta.median(mavalue)}')
            print()
        elif verbosity == 2:
            print(f'species: {j}')
            print(f'the mean of the species: {j} is: {sta.mean(mavalue)}')
            print(f'the mode of the species: {j} is: {sta.mode(mavalue)}')
            print(f'the median of the species: {j} is: {sta.median(mavalue)}')
            print()
    return print(f'the outliers of all pair are: {outliers}')

#df = pd.read_table(file_tsv)
#head = df.columns
#i, j = 0, 1
#oneline = df.iloc[i,j].strip()#.split(',') # [i,j] i = row; j = column


obj = read(file_tsv, count_dir)
out = statistics(obj)
print(out)




