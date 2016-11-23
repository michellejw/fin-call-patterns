# -*- coding: utf-8 -*-
"""
FIG_AX-KEMF.py

Script to generate figure showing Frequency vs IPI for only KEMF and Axial. 
This is the justification for extending the Axial data to an entire decade using
data from KENE

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
import pickle

df = pd.read_csv('../SEQ_CODE/ALL_seq_summary.csv')
sta = pd.read_csv('../SEQ_CODE/stationlist.csv')
# Force CI station names to match everything else
df.loc[df.station == 'J23A','station'] = 'CI-J23A'
df.loc[df.station == 'J06A','station'] = 'CI-J06A'
df.loc[df.station == 'J63A','station'] = 'CI-J63A'
df.loc[df.station == 'G03A','station'] = 'CI-G03A'
df.loc[df.station == 'G30A','station'] = 'CI-G30A'


df.datevec = pd.to_datetime(df.datevec)

df['year'] = df.datevec.dt.year
df['month'] = df.datevec.dt.month
unqsta = np.unique(df['station'])

# Set up colors (see cols.py)
coldict = pickle.load( open( "cols.p", "rb" ) )    
    
# loop through each row in df
stacols = []
for idx in range(len(df['station'])):
    stacols.append(coldict[df['station'][idx]])
    
df['stacolors'] = stacols


# Subset to obtain only KEMF, and Axial
df_decade = df[((df.station == 'KEMF') | 
            (df.station=='AX')) & (df.peakcounts > 50) 
            & (df.year > 2010) & (df.year < 2014)]
            

            
            
#----- PLOTTING -----*
plt.figure(1,figsize=(6,7))


# Subplot 1: IPI vs Frequency 
ax1 = plt.subplot2grid((4,1), (0,0), rowspan=3)
ax1.scatter(df_decade['freq'],df_decade['ipi'],s=df_decade['peakcounts']*.4,
        alpha=1,linewidths=2,
        edgecolors=df_decade['stacolors'],c='None')
ax1.set_xlabel('Frequency (Hz)')
ax1.set_ylabel('IPI (s)')

# dummy legend entries
axcol = plt.scatter([],[],s=100,alpha=1,linewidths=2,edgecolors=coldict['AX'],c='None',label='Axial')
#kenecol = plt.scatter([],[],s=100,alpha=1,linewidths=2,edgecolors=coldict['KENE'],c='None',label='KENE')
kemfcol = plt.scatter([],[],s=100,alpha=1,linewidths=2,edgecolors=coldict['KEMF'],c='None',label='KEMF')

pcounts = df_decade['propcounts']
pkcounts = df_decade['peakcounts']
s1 = plt.scatter([],[],s=100*.4,edgecolors='grey',c='None',alpha=2,label='100',linewidths=2)
s2 = plt.scatter([],[],s=300*.4,edgecolors='grey',c='None',alpha=2,label='300',linewidths=2)
s3 = plt.scatter([],[],s=800*.4,edgecolors='grey',c='None',alpha=2,label='800',linewidths=2)

legend1 = plt.legend(handles=[axcol,kemfcol],loc=3,bbox_to_anchor=(0.2,-0.35),scatterpoints=1,frameon=False,fontsize=11)
legend2 = plt.legend(handles=[s1,s2,s3],loc=3,bbox_to_anchor=(0.5,-0.4),scatterpoints=1,frameon=False,fontsize=11)
plt.gca().add_artist(legend1)


plt.savefig('final-figs/FIG_decade.tif')



            


