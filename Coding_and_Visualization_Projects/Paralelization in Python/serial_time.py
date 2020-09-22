#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 29 21:21:29 2020

@author: allsd
"""
import pandas as pd
import matplotlib.pyplot as plt
import time as t
import statistics
import re

pattern = input("Insert Pattern here (only characters): ") 

p = re.compile("[a-zA-Z]")

while(not bool(p.match(pattern))):
    pattern = input("Only characters from a to z accepted. Please repeat the pattern: ")

num = int(input("Insert number of iterations to generate a mean time: "))

pattern = pattern.upper()

def count_pattern(sequence, pattern):
    return(str(sequence).count(pattern))
    

time = []
proteins = pd.read_csv("proteins.csv", ";")

def serial(i, proteins, pattern) :
    
    start = t.time()
    
    data_dict = proteins.to_dict()
    
    data_dict['count']={}
    
    results = []
    
    for index, sequence in data_dict['sequence'].items(): 
        results.append(count_pattern(sequence, pattern))

    proteins['count'] = results
    
    end = t.time()
    time.append(end - start)
    
    print("Time of founding patterns in serial with dictionaries: " + str(end - start) + " seconds")
    
    return(proteins)
        
    
for i in range(num):
    proteins = serial(i, proteins, pattern)

 
occurences = proteins[proteins["count"] != 0]

plt.hist(occurences["count"],bins = 25)
    
print(occurences.loc[occurences['count'].idxmax()]['sequence'])
    
    
print("The mean time of serial is " + str(statistics.mean(time)) + " seconds")
