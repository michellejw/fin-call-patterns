# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 11:16:57 2016

@author: michw
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
#import statsmodels.api as sm
#from palettable.tableau import GreenOrange_12
import pickle

df = pd.read_csv('../SEQ_CODE/ALL_seq_summary.csv')
sta = pd.read_csv('../SEQ_CODE/stationlist.csv')

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


# Subset to obtain only KENE, KEMF, and Axial
df_decade = df[((df.station == 'KENE') | (df.station == 'KEMF') | 
            (df.station=='AX')) & (df.peakcounts > 50)]
            
#colors = {'KENE':'red','KEMF':'blue','AX':'orange'} 
#df_decade['stacols'] = df_decade['station'].apply(lambda x: colors[x])          
            
            
            
            
#----- PLOTTING -----*
plt.figure(1,figsize=(9,5))


# Subplot 1: IPI vs Frequency 
ax1 = plt.subplot2grid((4,2), (0,0), rowspan=3)
ax1.scatter(df_decade['freq'],df_decade['ipi'],s=df_decade['peakcounts']*.4,
        alpha=1,linewidths=2,
        edgecolors=df_decade['stacolors'],c='None')
ax1.set_xlabel('Frequency (Hz)')
ax1.set_ylabel('IPI (s)')

# dummy legend entries
axcol = plt.scatter([],[],s=100,alpha=1,linewidths=2,edgecolors=coldict['AX'],c='None',label='Axial')
kenecol = plt.scatter([],[],s=100,alpha=1,linewidths=2,edgecolors=coldict['KENE'],c='None',label='KENE')
kemfcol = plt.scatter([],[],s=100,alpha=1,linewidths=2,edgecolors=coldict['KEMF'],c='None',label='KEMF')

pcounts = df_decade['propcounts']
pkcounts = df_decade['peakcounts']
s1 = plt.scatter([],[],s=100*.4,edgecolors='grey',c='None',alpha=2,label='100',linewidths=2)
s2 = plt.scatter([],[],s=300*.4,edgecolors='grey',c='None',alpha=2,label='300',linewidths=2)
s3 = plt.scatter([],[],s=800*.4,edgecolors='grey',c='None',alpha=2,label='800',linewidths=2)

legend1 = plt.legend(handles=[axcol,kenecol,kemfcol],loc=3,bbox_to_anchor=(0.2,-0.43),scatterpoints=1,frameon=False,fontsize=11)
legend2 = plt.legend(handles=[s1,s2,s3],loc=3,bbox_to_anchor=(0.5,-0.43),scatterpoints=1,frameon=False,fontsize=11)
plt.gca().add_artist(legend1)

# Subplot 2: IPI vs year
# linear regression fit
df_main = df_decade[(df_decade['ipi'] > 21)&(df_decade['freq']<21)]
mod_ipi= smf.ols('ipi~year',data=df_main).fit()
ax2 = plt.subplot2grid((4,2), (0,1),rowspan=2)
ax2.yaxis.tick_right()
ax2.yaxis.set_label_position('right')
ax2.set_ylabel('IPI (s)',rotation=270,labelpad=15)
ax2.set_xticklabels([])
ax2.scatter(df_decade['year'],df_decade['ipi'],s=df_decade['peakcounts']*.4,
            alpha=1,linewidths=2,
            edgecolors=df_decade['stacolors'],c='None')
ax2.plot(np.array([2002,2014]),np.array([2002,2014])*mod_ipi.params[1] + mod_ipi.params[0],c='darkgray')
ax2.set_xlim([2002,2014])

# Subplot 3: Frequency vs year
mod_freq = smf.ols('freq~year',data=df_main).fit()
ax3 = plt.subplot2grid((4,2),(2,1),rowspan=2)
ax3.yaxis.tick_right()
ax3.yaxis.set_label_position('right')
ax3.set_ylabel('Frequency (Hz)',rotation=270,labelpad=15)
ax3.set_yticklabels(np.append(np.arange(17,25,1),''))
ax3.set_xticks(np.arange(2004,2014,2))
ax3.set_xlim([2002,2014])
ax3.set_xlabel('Year')
ax3.scatter(df_decade['year'],df_decade['freq'],s=df_decade['peakcounts']*.4,
            alpha=1,linewidths=2,
            edgecolors=df_decade['stacolors'],c='None')
ax3.plot(np.array([2002,2014]),np.array([2002,2014])*mod_freq.params[1] + mod_freq.params[0],c='darkgray')
ax3.set_xlim([2002,2014])


plt.tight_layout(h_pad=-2,w_pad=0.1)
#plt.tight_layout()
#plt.tight_layout(w_pad=0.1)

plt.savefig('FIG_decade.tif')



            