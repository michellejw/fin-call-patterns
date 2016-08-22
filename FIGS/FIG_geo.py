# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 15:20:50 2016

@author: michw
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle

df = pd.read_csv('../SEQ_CODE/ALL_seq_summary.csv')
sta = pd.read_csv('../SEQ_CODE/stationlist.csv')

df.datevec = pd.to_datetime(df.datevec)
df['year'] = df.datevec.dt.year
df['month'] = df.datevec.dt.month
df['dt'] = df.datevec.dt.date
unqsta = np.unique(df['station'])

# Set up colors (see cols.py)
coldict = pickle.load( open( "cols.p", "rb" ) )    

# loop through each row in df
stacols = []
for idx in range(len(df['station'])):
    stacols.append(coldict[df['station'][idx]])
    
df['stacolors'] = stacols

# Subset the data            
df = df[(df['datevec'] > '2011-07-01') & (df['datevec'] < '2012-06-30')]


# Subset for peakcounts
df = df[(df.peakcounts > 50)]

# ----- Plotting ----- #

f1 = plt.figure(2,figsize=(6,8))

ax1 = plt.subplot2grid((4,1), (0,0), rowspan=3)
ax1.scatter(df['freq'],df['ipi'],s=df['peakcounts'],
        linewidths=2,edgecolors=df['stacolors'],c='None')
ax1.set_xlabel('Frequency (Hz)')
ax1.set_ylabel('IPI (s)')

# dummy legend entries
legend_cols = []
legend_size = []
for idx in np.arange(len(unqsta)):
    legend_cols.append(plt.scatter([],[],s=100, linewidths=2,c='None',edgecolors=coldict[unqsta[idx]],label=unqsta[idx]))


pcounts = df['propcounts']
pkcounts = df['peakcounts']
s1 = plt.scatter([],[],s=100*.4,edgecolors='grey',c='None',linewidths=1,label='100')
s2 = plt.scatter([],[],s=300*.4,edgecolors='grey',c='None',linewidths=1,label='300')
s3 = plt.scatter([],[],s=800*.4,edgecolors='grey',c='None',linewidths=1,label='800')

legend1 = plt.legend(handles=legend_cols,title='Stations',ncol=3,loc=3,bbox_to_anchor=(-.05,-0.34),scatterpoints=1,frameon=False,fontsize=11)
#legend1 = plt.legend(handles=[axcol,kenecol,kemfcol],loc=3,bbox_to_anchor=(0.2,-0.43),scatterpoints=1,frameon=False,fontsize=11)
legend2 = plt.legend(handles=[s1,s2,s3],title='Counts',loc=3,bbox_to_anchor=(0.8,-0.34),scatterpoints=1,frameon=False,fontsize=11)
plt.gca().add_artist(legend1)

plt.savefig('FIG_geo.tif')
