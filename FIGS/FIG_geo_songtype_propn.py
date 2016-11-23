# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 14:54:14 2016

@author: michw
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
label_size = 14
mpl.rcParams['xtick.labelsize'] = label_size 
mpl.rcParams['ytick.labelsize'] = label_size
mpl.rcParams['axes.labelsize'] = label_size # might only need this one?


df = pd.read_csv('../SEQ_CODE/ALL_seq_summary4.csv')
sta = pd.read_csv('../SEQ_CODE/stationlist.csv')
sta = sta.sort_values(by='latitude',ascending=False)
df = pd.merge(df,sta,on='station')

df.datevec = pd.to_datetime(df.datevec)
df['year'] = df.datevec.dt.year
df['month'] = df.datevec.dt.month
df['dt'] = df.datevec.dt.date
unqsta = np.unique(df['station'])
# create month-year variable
df['monthyear'] = df['datevec'].apply(lambda x: x.strftime('%b-%Y')) 

df['songtype'] = df.ipi>=21 # True=singlet, False=doublet
df['notetype'] = df.freq>=22 # True=Note B, False=Note A

# Set up plot region       
f1, axarr = plt.subplots(5,1, figsize=(12,14))
plt.setp([a.get_xticklabels() for a in axarr[:-1]], visible=False)
f1.subplots_adjust(hspace=0.15,wspace=0.05)

monthlist = np.array(('Nov-2011','Dec-2011','Jan-2012','Feb-2012','Mar-2012'))


for m in np.arange(len(monthlist)):
    dfsub = df[(df.peakcounts > 50) & (df.monthyear==monthlist[m])]
    # groupby station and song type (doublet/singlet)
    G_sta_song = dfsub.groupby(['station','songtype']).peakcounts.sum()
    G_sta = dfsub.groupby(['station']).peakcounts.sum()
    G_pct = G_sta_song.div(G_sta,level='station')
    
    G_pct2 = G_pct.unstack(level=1)
    # list of stations in desired order (excluding kene and cz)
    stalist = list(sta.station[(sta.station != 'KENE') & (sta.station != 'CZ')])
    G_pct2 = G_pct2.reindex(stalist)
    # Nans are actually places where there were no songs of a particular type
    G_pct2 = G_pct2.fillna(0)
    
    # plot bar chart
    G_pct2[0].plot.barh(ax=axarr[m],color='royalblue') # doublets
    G_pct2[1].plot.barh(left=G_pct2[0],color='orange',ax=axarr[m]) # singlet
    
    axarr[m].set_ylabel(monthlist[m])

axarr[m].set_xlabel('Proportion of counts')

# dummy legend entries
doublet = plt.scatter([], [],marker='o',s=100,alpha=1,color="royalblue",label='Doublet');
singlet = plt.scatter([], [],marker='o',s=100,alpha=1,color="orange",label='Singlet');

plt.legend(handles=[doublet,singlet],loc=3,bbox_to_anchor=(-0.04,-0.6),fontsize=14)

plt.savefig('final-figs/FIG_geo_songtype_propn.png')
