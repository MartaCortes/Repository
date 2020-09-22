#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 29 21:40:01 2020

@author: allsd
"""

import pandas as pd
import matplotlib.pyplot as plt
import multiprocessing as mp
import statistics
import time as t
import re

p = re.compile("[a-zA-Z]")

pattern = input("Insert Pattern here (only characters): ") 

while(not bool(p.match(pattern))):
    pattern = input("Only characters from a to z accepted. Please repeat the pattern: ")

num = int(input("Insert number of iterations to generate a mean time: "))

pattern = pattern.upper()
proteins = pd.read_csv("proteins.csv", ";")

def count_pattern(sequence, pattern):
        return(str(sequence).count(pattern))
        
time = []

def multi (i, proteins, pattern):
    start = t.time()
    
    data_dict = proteins.to_dict()
    
    pool = mp.Pool(mp.cpu_count())
    
    array_of_args = [(sequence, pattern) for sequence in data_dict['sequence'].values()]
    
    results = pool.starmap(count_pattern, array_of_args)
      
    pool.close()
    
    proteins['count'] = results
        
    end = t.time()
    
    print("Time of founding patterns in multiprocessing with dictionaries: " + str(end - start) + " seconds")
    
    time.append(end - start)
    
    return(proteins)

for i in range(num):
    proteins = multi(i, proteins, pattern)

 
occurences = proteins[proteins["count"] != 0]

plt.hist(occurences["count"],bins = 25)
    
print(occurences.loc[occurences['count'].idxmax()]['sequence'])
    
    
print("The mean time of multiprocessing is " + str(statistics.mean(time)) + " seconds")