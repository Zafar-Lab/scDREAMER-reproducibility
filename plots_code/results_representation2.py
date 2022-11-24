
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Created on Tue Jan 11 12:37:53 2022
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


print("Current Working Directory " , os.getcwd())
os.chdir("/Users/ajitashree/Documents/OneDrive - IIT Kanpur/data-integration/")
print("Now Current Working Directory " , os.getcwd())

df = pd.read_excel('./All_metrics_15_Sep.xlsx', sheet_name = 'all_metrics_revision') #scDREAMER++, plot

#name = 'Human mouse'
#name = 'Immune human'
#name = 'Human mouse'
#name = "Lung"
name = "Immune human"

'''
name = 'Pancreas'  - final
name = 'Lung' - v final
name = 'Immune human' - v final
name = 'Human mouse' #** - v final
name = 'Human retina' - final
'''

"""
scVI - cyan
Harmony - orange
Seurat - magenta
BBKNN - purple
Scanorama - yellow
INSCT - blue
iMAP - brown
Liger - light green
fastMNN - light pink
scANVI - removed from main figure
scDREAMER - Green
scDREAMER++ - Red

"""
#df_ = hum_mou # sim1, pan, lung, sim1, retina, imm_hum_mou, retina, tabula, hum_mou , imm_hum
df_ = df[df['Dataset'] == name].reset_index(drop = True)



dct = {'scVI' : '#28DDED', 'Seurat' : '#994363', 'Harmony': '#ED7A28',
       'BBKNN': '#B626D3', 'Scanorama': '#EDBF28', 'scDREAMER': '#086E28',
       'INSCT' : '#286CED'}


dct = {'scVI' : '#28DDED', 'Seurat' : '#994363', 'Harmony': '#ED7A28',
       'BBKNN': '#B626D3', 'Scanorama': '#EDBF28', 
       'INSCT' : '#286CED', 'fastMNN':  "#FFB6C1", "iMAP" : "#964B00",
       'LIGER' : '#90EE90',
       'scDREAMER': '#086E28'
       #,'scDREAMER++' : "#ff0000"
       }

"""
scDREAMER, scDREAMER-woDis
scDREAMER-woBC

dct = {"scDREAMER": "#086E28",
"scDREAMER-woDis" : "#994363",
"scDREAMER-woBC" : "#28DDED",
"scDREAMER++ V2 - crossE" : '#ED7A28',
"scANVI" : '#B626D3',
"scDREAMER++ w/o BC or Dis":'#EDBF28',
"scDREAMER++ w/o BC or Dis V2" : "#0000FF",
"scgen" : "#964B00",
"scDREAMER++ V2 - crossEV2" : "#ff0000"
}

"""
lst_mthds = dct.keys()
df_ = df_[df_['Method'].isin(lst_mthds)]

clr = []
for index, row in df_.iterrows():
    clr += [dct[row['Method']]]
    
df_ = df_[df_.columns[:14]]
df_ = df_.fillna(0)

# below steps for other kind of datasets
#new_header = df_.iloc[0] #grab the first row for the header
#df_ = df_[1:] #take the data less the header row
#df_.columns = new_header



### cLISI vs iLISI plots ########
#grid()

colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', 
                   '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

df_['color'] = clr #colors[:len(df_)]

'''
for index, row in df_.iterrows():
    
    #print (row)
    plt.plot(float(row['iLISI']), float(row['cLISI']), 'bo', marker = '.', ms = 14, color = row['color'], label = row['Method'])


mi = df_['cLISI'].min()
mx = df_['cLISI'].max()

ni = df_['iLISI'].min()
nx = df_['iLISI'].max()

ylim(mi - 0.01, mx + 0.01) 
xlim(ni - 0.05, nx + 0.05) 

plt.legend(loc='center left', bbox_to_anchor = (1, 0.5)) 
#bbox_to_anchor=(0.5, 1.05); (1.1, 1.05)
#plt.title(name + ' cLISI vs iLISI', fontsize = 15)
plt.xlabel('iLISI', fontsize = 10)
plt.ylabel('cLISI', fontsize = 10)
show()
'''
#isolated sc

def plot_bar(df_, col_name):
    
    rc('axes', linewidth = 2)

    if (name == 'Human mouse' and col_name == 'kBET'):
        df_ = df_[df_['Method'] != 'BBKNN'].reset_index(drop = True)
        
    if (col_name == 'ASW label' or col_name == 'ASW label/batch' or col_name == 'PCR batch' or col_name == 'isolated silhouette coefficient'):
        df_ = df_[df_['Method'] != 'BBKNN'].reset_index(drop = True)



    fig = plt.figure(figsize = (6, 4))
    ax = df_[col_name].plot(kind="bar", color = df_['color'])

    barWidth = 0.3
    br1 = np.arange(len(df_))
    
    
    """
    plt.figure(figsize=(12, 8))
    ax = freq_series.plot(kind="bar")
    ax.set_title("Amount Frequency")
    ax.set_xlabel("Amount ($)")
    ax.set_ylabel("Frequency")
    ax.set_xticklabels(x_labels)
    
    rects = ax.patches
    
    # Make some labels.
    labels = [f"label{i}" for i in range(len(rects))]
    
    for rect, label in zip(rects, labels):
        height = rect.get_height()
        ax.text(
            rect.get_x() + rect.get_width() / 2, height + 5, label, ha="center", va="bottom"
        )
    
    plt.show()
    
    """
    rects = ax.patches
    ax.set_xticklabels(df_['Method'], rotation = 60, fontname='Arial', fontsize = 10)
    
    #barWidth = 0.1
    #ax.bar(br1, df_[col_name], color = df_['color'], width = barWidth,
    #         label = col_name)
        
    for rect, label in zip(rects, df_[col_name]):
        height = rect.get_height()
        ax.text(
        rect.get_x() + rect.get_width() / 2, height, round(label, 2), ha="center", va="bottom"
        )
    

    
    mi = df_[col_name].min()
    mx = df_[col_name].max()
    
    ylim(mi - 0.01, min(mx*1.05, 1.0)) #(0.6, 1)
    ylim(mi - 0.01, min(mx*1.05, 1.0))
    ax.legend()
    
    #ax = gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(14)
        tick.label1.set_fontweight('bold')
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(14)
        tick.label1.set_fontweight('bold')
        
        
    plt.subplot(111).spines['right'].set_visible(False)
    plt.subplot(111).spines['top'].set_visible(False)
        
    plt.ylabel(col_name, fontsize = 15, fontname='Arial', fontweight = 'bold')
    #plt.xlabel("Methods", fontsize = 15, fontname='Arial', fontweight = 'bold')
    #plt.title(name + ": " + col_name, fontsize = 20, fontweight = 'bold')
    
    plt.tight_layout()
    plt.gcf().subplots_adjust(bottom = 0.3)
    #plt.gcf().subplots_adjust(top = 3)
    #plt.gcf().subplots_adjust(top=5)

    """
    path = "/Users/ajitashree/Desktop/Projects/Project_DI/_final_plots_mar22/Figures/"
    
    if "/" not in col_name:
        plt.savefig(path + name + '_'+ col_name + '.png')
    else:
        plt.savefig(path + name + '_'+ col_name.split(' ')[0] + '.png', bbox_inches = "tight")
    """
    plt.show()
    
def plot_bar_(df_, col_name):
    
    rc('axes', linewidth = 2)

    if (name == 'Human mouse' and col_name == 'kBET'):
        df_ = df_[df_['Method'] != 'BBKNN'].reset_index(drop = True)
        
    if (col_name == 'ASW label' or col_name == 'ASW label/batch' or col_name == 'PCR batch' or col_name == 'isolated silhouette coefficient'):
        df_ = df_[df_['Method'] != 'BBKNN'].reset_index(drop = True)


    fig = plt.figure(figsize = (6, 4))
    #grid()
    barWidth = 0.3
    br1 = np.arange(len(df_))
    
    plt.bar(br1, df_[col_name], color = df_['color'], width = barWidth,
             label = col_name)
    
    plt.xticks([r for r in range(len(df_))],
            df_['Method'], rotation = 90, fontname='Arial', fontsize = 10)
    
    #xlim(0.6, 1)
    mi = df_[col_name].min()
    mx = df_[col_name].max()
    
    ylim(mi - 0.01, min(mx*1.05, 1.0)) #(0.6, 1)
    ylim(mi - 0.01, min(mx*1.05, 1.0))
    #plt.legend()
    
    ax = gca()
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(14)
        tick.label1.set_fontweight('bold')
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(14)
        tick.label1.set_fontweight('bold')
        
        
    plt.subplot(111).spines['right'].set_visible(False)
    plt.subplot(111).spines['top'].set_visible(False)
        
    plt.ylabel(col_name, fontsize = 15, fontname='Arial', fontweight = 'bold')
    #plt.xlabel("Methods", fontsize = 15, fontname='Arial', fontweight = 'bold')
    #plt.title(name + ": " + col_name, fontsize = 20, fontweight = 'bold')
    
    plt.tight_layout()
    plt.gcf().subplots_adjust(bottom = 0.3)
    #plt.gcf().subplots_adjust(top = 3)
    #plt.gcf().subplots_adjust(top=5)

    """
    path = "/Users/ajitashree/Desktop/Projects/Project_DI/_final_plots_mar22/Figures/"
    
    if "/" not in col_name:
        plt.savefig(path + name + '_'+ col_name + '.png')
    else:
        plt.savefig(path + name + '_'+ col_name.split(' ')[0] + '.png', bbox_inches = "tight")
    """
    plt.show()
    
#### isolated f1

for i in df_.columns[2:-1]:
    print (i)
    plot_bar(df_, i)
  


### batch correction ############
'''
fig = plt.figure(figsize = (10, 5))
#grid()
# creating the bar plot
##plt.bar(df_['NMI_cluster/label'], df_['ARI_cluster/label'],width = 0.4)


barWidth = 0.15
br1 = np.arange(len(df_))
br2 = [x + barWidth for x in br1]
br3 = [x + barWidth for x in br2]
br4 = [x + barWidth for x in br3]

plt.bar(br1, df_['ASW label/batch'], color ='r', width = barWidth,
        edgecolor ='grey', label ='ASW')

plt.bar(br2, df_['PCR batch'], color ='b', width = barWidth,
        edgecolor ='grey', label ='PCR batch')

plt.bar(br3, df_['graph connectivity'], color ='y', width = barWidth,
        edgecolor ='grey', label ='graph connectivity')

plt.bar(br4, df_['kBET'], color ='g', width = barWidth,
        edgecolor ='grey', label ='kbet')

plt.xticks([r + barWidth for r in range(len(df_))],
        df_['Method'], rotation = 60)

plt.legend()
#plt.legend(loc='center left', bbox_to_anchor = (1, 0.5))
plt.ylabel("NMI, ARI")
plt.xlabel("Methods")
#plt.title(name + ": batch correction metrics", fontsize = 15)

plt.show()
'''

data = df_.copy(deep = True)
# Scaling method 1: 
d = data.drop(['Dataset', 'Method', 'color'], axis=1)
#df_ = (d-d.min())/(d.max()-d.min())
#df_ = pd.concat((df_, data[['Dataset', 'Method', 'color']]), 1)
 
#print("Scaled Dataset Using Pandas")
#df_.head()

####################

'''
fig = plt.figure(figsize = (10, 5))
#grid()
# creating the bar plot
##plt.bar(df_['NMI_cluster/label'], df_['ARI_cluster/label'],width = 0.4)


barWidth = 0.15
br1 = np.arange(len(df_))
br2 = [x + barWidth for x in br1]
br3 = [x + barWidth for x in br2]
br4 = [x + barWidth for x in br3]

plt.bar(br1, df_['ASW label/batch'], color ='r', width = barWidth,
        edgecolor ='grey', label ='ASW')

plt.bar(br2, df_['PCR batch'], color ='b', width = barWidth,
        edgecolor ='grey', label ='PCR_batch')

plt.bar(br3, df_['graph connectivity'], color ='y', width = barWidth,
        edgecolor ='grey', label ='graph connectivity')

plt.bar(br4, df_['kBET'], color ='g', width = barWidth,
        edgecolor ='grey', label ='kbet')

plt.xticks([r + barWidth for r in range(len(df_))],
        df_['Method'], rotation = 45)

plt.legend()
#plt.legend(loc='center left', bbox_to_anchor = (1, 0.5))
plt.ylabel("NMI, ARI")
plt.xlabel("Methods")
#plt.title(name + ": batch correction metrics", fontsize = 15)

plt.savefig('test.pdf')
plt.show()
'''


### using min max scaler from sklearn ######

from sklearn.preprocessing import MinMaxScaler
scaler = MinMaxScaler()

bbknn_name = "bbknn" # "BBKNN"

data = df_.copy(deep = True)
data.loc[data.Method == bbknn_name,'ASW label/batch'] = np.nan
data.loc[data.Method == bbknn_name,'ASW label'] = np.nan
data.loc[data.Method == bbknn_name,'PCR batch'] = np.nan
data.loc[data.Method == bbknn_name,'isolated silhouette coefficient'] = np.nan
d = data.drop(['Dataset', 'Method', 'color'], axis = 1)


M = df_.copy(deep = True)
df_ = scaler.fit_transform(d.to_numpy())
df_ = pd.DataFrame(df_, columns= d.columns)

data = data.reset_index(drop = True)

df_ = pd.concat((df_, data[['Dataset', 'Method', 'color']]), 1)

T = df_.copy(deep = True)

### Conclusion: both the scaling methods are giving same
# https://www.geeksforgeeks.org/how-to-scale-pandas-dataframe-columns/

#df_['Composite batch-correction score'] = df_[['ASW label/batch', 'PCR batch', 'graph connectivity', 'kBET']].mean(axis = 1)
#df_['bio conservation metrics average'] = df_[['NMI_cluster/label', 'ARI_cluster/label', 'ASW_label', 'isolated f1', 'isolated sc']].mean(axis = 1)
df_['Composite bio-conservation score'] = df_[['NMI cluster/label', 'ARI cluster/label', 'ASW label']].mean(axis = 1)
df_['Composite batch-correction score'] = df_[['ASW label/batch', 'PCR batch', 'graph connectivity', 'kBET']].mean(axis = 1)
df_['Composite isolated label score'] = df_[['isolated silhouette coefficient', 'isolated f1 score']].mean(axis = 1)


bbknn = df_[df_['Method'] == bbknn_name]

if (name != 'Human mouse'):
    avg = (bbknn['kBET']+ bbknn['graph connectivity'])/2
else:
    avg = (bbknn['graph connectivity'])

avg_ = (bbknn['NMI cluster/label']+ bbknn['ARI cluster/label'])/2
avg_i = bbknn['isolated f1 score']

df_.loc[df_.Method==bbknn_name,'Composite batch-correction score'] = avg
df_.loc[df_.Method==bbknn_name,'Composite bio-conservation score'] = avg_
df_.loc[df_.Method==bbknn_name,'isolated f1 score'] = avg_i

print (df_)

'''
if (name == 'Human mouse'):
    bbknn = df_[df_['Method'] == 'BBKNN']
    avg = (bbknn['ASW label/batch'] + bbknn['PCR batch'] + bbknn['graph connectivity'])/3
    df_.loc[df_.Method=='BBKNN','Composite batch-correction score'] = avg
'''
  

df_['Combined composite score'] = df_[['Composite bio-conservation score', 'Composite batch-correction score']].mean(axis = 1)


plot_bar(df_, 'Composite batch-correction score')
plot_bar(df_, 'Composite bio-conservation score')
plot_bar(df_, 'Composite isolated label score')
plot_bar(df_, 'Combined composite score')


val = df_['Combined composite score'][-1:]
df_['combined comp % improvement'] = (100* (val.values[0] - df_['Combined composite score'])/df_['Combined composite score'])

val = df_['Composite bio-conservation score'][-1:]
df_['bio comp % improvement'] = (100* (val.values[0] - df_['Composite bio-conservation score'])/df_['Composite bio-conservation score'])

val = df_['Composite batch-correction score'][-1:]
df_['batch comp % improvement'] = (100* (val.values[0] - df_['Composite batch-correction score'])/df_['Composite batch-correction score'])

val = df_['Composite isolated label score'][-1:]
df_['iso comp % improvement'] = (100* (val.values[0] - df_['Composite isolated label score'])/df_['Composite isolated label score'])

writer = pd.ExcelWriter('./' + name + 'ablation_comparison.xlsx', engine = 'xlsxwriter')
df_.to_excel(writer, sheet_name = name)
writer.save()
writer.close()



##############################################################################
import scanpy as sc

#/Users/ajitashree/Downloads
#path = './Downloads/' 
#ann1 = sc.read_h5ad(path + 'HCL_Fig1_adata.h5ad')
#ann2 = sc.read_h5ad(path + 'MCA_BatchRemoved_Merge_dge.h5ad.h5ad')
























