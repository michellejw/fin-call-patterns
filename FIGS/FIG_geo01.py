# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 14:54:14 2016

@author: michw
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle

df = pd.read_csv('../SEQ_CODE/ALL_seq_summary4.csv')
sta = pd.read_csv('../SEQ_CODE/stationlist.csv')
sta = sta.sort_values(by='latitude',ascending=False)
df = pd.merge(df,sta,on='station')

df.datevec = pd.to_datetime(df.datevec)
df['year'] = df.datevec.dt.year
df['month'] = df.datevec.dt.month
df['dt'] = df.datevec.dt.date
unqsta = np.unique(df['station'])

df['songtype'] = df.ipi>=21
df['notetype'] = df.freq>=22

# Set up colors (see cols_months.py)
coldict = pickle.load( open( "cols_months.p", "rb" ) )    

# loop through each row in df
monthcols = []
for idx in range(len(df['month'])):
    monthcols.append(coldict[df['month'][idx]])
    
df['monthcols'] = monthcols

# Subset the data by date and peakcounts
dfALL = df[(df.peakcounts >50)]
df = df[(df.peakcounts > 50) & (df.datevec >= '2011-11-01') & (df.datevec <= '2012-03-01')]
df = df.sort_values(by='latitude')

# groupby station and song type (doublet/singlet)
G_sta_song = df.groupby(['station','songtype']).songtype.count()
G_sta = df.groupby(['station']).songtype.count()
G_pct = G_sta_song.div(G_sta,level='station')

G_pct2 = G_pct.unstack(level=1)
# list of stations in desired order (excluding kene and cz)
stalist = list(sta.station[(sta.station != 'KENE') & (sta.station != 'CZ')])
G_pct2 = G_pct2.reindex(stalist)

# plot bar chart
G_pct2[0].plot.barh()
G_pct2[1].plot.barh(left=G_pct2[0],color='magenta')


# More plotting
# groupby station and song type (doublet/singlet)
G_sta_song = df.groupby(['station','songtype']).peakcounts.sum()
G_sta = df.groupby(['station']).peakcounts.sum()
G_pct = G_sta_song.div(G_sta,level='station')

G_pct2 = G_pct.unstack(level=1)
# list of stations in desired order (excluding kene and cz)
stalist = list(sta.station[(sta.station != 'KENE') & (sta.station != 'CZ')])
G_pct2 = G_pct2.reindex(stalist)

# plot bar chart
plt.figure(33)
G_pct2[0].plot.barh()
G_pct2[1].plot.barh(left=G_pct2[0],color='orange')


# plot ipi vs latitude
dfALLsub = dfALL[dfALL.freq < 20]
dfsub = df[df.freq < 20]
plt.figure(1)
plt.clf()
#plt.scatter(dfALLsub.ipi,dfALLsub.latitude,s=dfALLsub.peakcounts,alpha=0.2,c='lightgrey')
plt.scatter(dfsub.ipi,dfsub.latitude)
plt.grid()

