# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np

with open("Data_7505.txt") as inp:
    rw = list(zip(*(line.strip().split('\t') for line in inp)))
    data = list(rw[:][:]);
with open("Reporters_7505.txt") as inp:
    rep = list(zip(*(line.strip('\n').split('\t') for line in inp)))

i = rep[5].index('NR_015450.1')
N = rep[4][i]
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