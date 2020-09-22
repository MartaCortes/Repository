# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas as pd
import matplotlib.pyplot as plt

import multiprocessing as mp

import time as t

print("Number of processors: ", mp.cpu_count())

pattern = input("Insert Pattern here (only characters): ") 

pattern = pattern.upper()

proteins = pd.read_csv("proteins.csv", ";")

proteins["count"] = 0

start = t.time()

pool = mp.Pool(mp.cpu_count())

def count_pattern(sequence, pattern):
    return (sequence.count(pattern))

proteins["count"] = [pool.apply(count_pattern, args=(str(proteins.iloc[i][1]),pattern)) for i in range(len(proteins.index))] 

pool.close()

#proteins["count"]= pool.map(count_pattern, [str(proteins.iloc[i][1]) for i in range(len(proteins.index))])

end = t.time()

print(end - start)

print("Time of founding patterns: " + str(end - start) + "seconds")

occurences = proteins[proteins["count"] != 0]

plt.hist(occurences["count"],bins = 25)

print(occurences.loc[occurences['count'].idxmax()]['sequence'])

