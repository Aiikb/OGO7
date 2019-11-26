# -*- coding: utf-8 -*-
"""
Use https://www.ncbi.nlm.nih.gov/nuccore/ to find nucleotide accession
"""
import matplotlib.pyplot as plt
import numpy as np

with open("Data_7177.txt") as inp:
    rw = list(zip(*(line.strip().split('\t') for line in inp)))
    data = list(rw[:][:]);
with open("Reporters_7177.txt") as inp:
    rep = list(zip(*(line.strip('\n').split('\t') for line in inp)))

i = rep[1].index('NM_198317')
N = rep[0][i]
j = data[0].index(N)    

for i in range(0,5):  
    print(data[i][0:5]) 

y = [] 
x = [] 
for i in range(1,len(data)):  
    print(data[i][j]) 
    y.append(data[i][j])
    x.append(data[i][0])
xp = np.arange(len(x))
plt.bar(xp,y)
plt.show()

