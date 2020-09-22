#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 18:17:36 2020

@author: allsd
"""

import pandas as pd
import matplotlib.pyplot as plt

import threading

import time as t

pattern = input("Insert Pattern here (only characters): ") 

pattern = pattern.upper()

proteins = pd.read_csv("proteins.csv", ";")

data_dict = proteins.to_dict()

def count_pattern(data_dict,pattern):
    data_dict['count']={}
    for key in data_dict['sequence']: 
        data_dict['count'][key] = str(data_dict['sequence'][key]).count(pattern)
        
#def count_pattern(proteins,pattern):
#    for i in range(len(proteins.index)): 
#        proteins.iloc[i,2] = str(proteins.iloc[i][1]).count(pattern)
# 
#def count_pattern(proteins,pattern):
#    for i in range(len(proteins.index)):
#        id_pt = proteins.iloc[i,0]
#        count_pt = str(proteins.iloc[i][1]).count(pattern)
#        data_dict[id_pt] = count_pt
#    return data_dict

start = t.time()

NUM_THREADS = 3

for num_thread in range(NUM_THREADS):
    th = threading.Thread(name='th%s' %num_thread, 
                         target= count_pattern,
                         args=(data_dict,pattern,))
    th.start()

for num_thread in range(NUM_THREADS):
    th.join()

end = t.time()

print(end - start)

result = pd.DataFrame.from_dict(data_dict)

print("Time of founding patterns: " + str(end - start) + " seconds")

occurences = result[result["count"] != 0]

plt.hist(occurences["count"],bins = 25)

print(occurences.loc[occurences['count'].idxmax()]['sequence'])
