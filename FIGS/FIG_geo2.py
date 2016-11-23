# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 15:20:50 2016

@author: michw
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle

#df = pd.read_csv('../SEQ_CODE/ALL_seq_summary.csv')
df = pd.read_csv('../SEQ_CODE/ALL_seq_summary4.csv')
sta = pd.read_csv('../SEQ_CODE/stationlist.csv')

df.datevec = pd.to_datetime(df.datevec)
df['year'] = df.datevec.dt.year
df['month'] = df.datevec.dt.month
df['dt'] = df.datevec.dt.date

# Subset the data            
df = df[(df['datevec'] > '2007-07-01') & (df['datevec'] < '2009-06-30')]

# Set up colors (see cols.py)
coldict = pickle.load( open( "cols.p", "rb" ) )    

# loop through each row in df
stacols = []
for idx in df.index:
    stacols.append(coldict[df['station'][idx]])
    
df['stacolors'] = stacols




# Subset for peakcounts
#df = df[(df.peakcounts > 50)]

# Unique stations after filtering
unqsta = np.unique(df['station'])

# ----- Plotting ----- #

f1 = plt.figure(5,figsize=(8,6))

ax1 = plt.subplot2grid((1,5), (0,0), colspan=4)
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

legend1 = plt.legend(handles=legend_cols,title='Stations',bbox_to_anchor=(1.25,0.7),scatterpoints=1,frameon=False,fontsize=11)
#legend1 = plt.legend(handles=[axcol,kenecol,kemfcol],loc=3,bbox_to_anchor=(0.2,-0.43),scatterpoints=1,frameon=False,fontsize=11)
legend2 = plt.legend(handles=[s1,s2,s3],title='Counts',bbox_to_anchor=(1.25,0.5),scatterpoints=1,frameon=False,fontsize=11)
plt.gca().add_artist(legend1)

plt.savefig('FIG_geo_axcz.tif')



# percentage of doublet calls at Colza
df_CZS_1 = df[(df['station'] == 'CZ') & (df['ipi']>=22) & (df['datevec'] <= '2008-03-30')]
df_CZD_1 = df[(df['station'] == 'CZ') & (df['ipi']<22) & (df['datevec'] <= '2008-03-30')]
y1d = np.sum(df_CZD_1['peakcounts'])
y1s = np.sum(df_CZS_1['peakcounts'])

df_CZS_2 = df[(df['station'] == 'CZ') & (df['ipi']>=22) & (df['datevec'] >= '2008-11-01')]
df_CZD_2 = df[(df['station'] == 'CZ') & (df['ipi']<22) & (df['datevec'] <= '2008-11-01')]
y2d = np.sum(df_CZD_2['peakcounts'])
y2s = np.sum(df_CZS_2['peakcounts'])

df_CZS_1 = df[(df['station'] == 'CZ') & (df['ipi']>=22) & (df['datevec'] <= '2008-03-30')]
df_CZD_1 = df[(df['station'] == 'CZ') & (df['ipi']<22) & (df['datevec'] <= '2008-03-30')]
y1d = np.sum(df_CZD_1['peakcounts'])
y1s = np.sum(df_CZS_1['peakcounts'])

df_CZS_2 = df[(df['station'] == 'CZ') & (df['ipi']>=22) & (df['datevec'] >= '2008-11-01')]
df_CZD_2 = df[(df['station'] == 'CZ') & (df['ipi']<22) & (df['datevec'] <= '2008-11-01')]
y2d = np.sum(df_CZD_2['peakcounts'])
y2s = np.sum(df_CZS_2['peakcounts'])


# doublet frequency
y1doublet_f = df_CZD_1[(df_CZD_1['freq']>20)&(df_CZD_1['freq']<22)]
meanCZdoubletfreq = np.mean(y1doublet_f['freq'])

y1doublet_f = df_AXD_1[(df_CZD_1['freq']>20)&(df_CZD_1['freq']<22)]
