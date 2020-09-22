# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas as pd
import matplotlib.pyplot as plt

pattern = input("Insert Pattern here (only characters): ") 

pattern = pattern.upper()

proteins = pd.read_csv("proteins.csv", ";")

proteins["count"] = 0

start = t.time()

def count_pattern(sequence, pattern):
    return (sequence.count(pattern))

for i in range(len(proteins.index)): 
    proteins.iloc[i,2] = count_pattern(str(proteins.iloc[i][1]), pattern)
 
end = t.time()

print("Time of founding patterns: " + str(end - start) + " seconds")
    
occurences = proteins[proteins["count"] != 0]

plt.hist(occurences["count"],bins = 25)

print(occurences.loc[occurences['count'].idxmax()]['sequence'])
