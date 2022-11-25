#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Thu Nov 24 16:27:34 2022
@author: ajitashree
"""

import pandas as pd
import numpy as np
import math
import pylab
from pylab import *
from sklearn import linear_model
from sklearn.metrics import mean_squared_error, r2_score
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 300

import os

from sklearn.preprocessing import MinMaxScaler
scaler = MinMaxScaler()


print("Current Working Directory " , os.getcwd())
os.chdir("/Users/ajitashree/Documents/OneDrive - IIT Kanpur/data-integration/")
print("Now Current Working Directory " , os.getcwd())

df = pd.read_excel('./scDREAMER++.xlsx') 
df_min = pd.read_excel('./scDREAMER++.xlsx', sheet_name = 'min')

name = "Lung" # Lung

df_ = df[df['Dataset'] == name].reset_index(drop = True)

dct = {'scgen' : '#28DDED', 'scanvi' : '#994363', 'scDREAMER++': '#086E28'}


data = df_.copy(deep = True)
d = data.drop(['Dataset', 'Method'], axis = 1)

#df_ = scaler.fit_transform(d.to_numpy())
#df_ = pd.DataFrame(df_, columns= d.columns)
#data = data.reset_index(drop = True)


I_min = d[d.columns].min()
I_max = d[d.columns].max()

T_min = df_min[df_min.columns].min() #d[d.columns].min()
T_max = I_min.copy(deep = True)

for c in d.columns[1:]:
    T_max[c] = 1

for c in d.columns[1:]:
  print (c)
  Numerator = (d[c] - I_min[c])*T_max[c] + T_min[c]*(I_max[c] - d[c])
  df_[c] = Numerator /(I_max[c] - I_min[c])
  
  
df_ = pd.concat((df_, data[['Dataset', 'Method']]), 1)


"""
I_min = out[out.columns].min()
I_max = out[out.columns].max()

T_min = inp[inp.columns].min()
T_max = inp[inp.columns].max()

out_scaled = out.copy(deep = True)
for c in out.columns:
  print (c)
  Numerator = (out[c] - I_min[c])*T_max[c] + T_min[c]*(I_max[c] - out[c])
  out_scaled[c] = Numerator /(I_max[c] - I_min[c])
  
"""
df_['Composite bio-conservation score'] = df_[['NMI cluster/label', 'ARI cluster/label', 'ASW label']].mean(axis = 1)
df_['Composite batch-correction score'] = df_[['ASW label/batch', 'PCR batch', 'graph connectivity', 'kBET']].mean(axis = 1)
df_['Composite isolated label score'] = df_[['isolated silhouette coefficient', 'isolated f1 score']].mean(axis = 1)
df_['Combined composite score'] = df_[['Composite bio-conservation score', 'Composite batch-correction score']].mean(axis = 1)



writer = pd.ExcelWriter('./New' + name + 'scdreamer_result__v2.xlsx', engine = 'xlsxwriter')
df_.to_excel(writer, sheet_name = name)
writer.save()
writer.close()



    







