# -*- coding: utf-8 -*-
"""
FIG_decade2.py

Script to generate frequency vs. time and IPI vs time plots to show decadal trend

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


# Subset to obtain only KENE, KEMF, and Axial
df_decade = df[((df.station == 'KENE') | (df.station == 'KEMF') | 
            (df.station=='AX')) & (df.peakcounts > 50)]
            
df_A = df_decade[(df.freq>=21)]     
df_B = df_decade[(df.freq < 21)] 
            
            
            
#----- PLOTTING -----*
plt.figure(1,figsize=(8,6))
plt.clf()

# Subplot 2: IPI vs year
# linear regression fit
df_main = df_B[(df_B['ipi'] > 20)]

mod_ipi= smf.ols('ipi~year',data=df_main).fit()
ax2 = plt.subplot2grid((4,3), (0,0),rowspan=2,colspan=2)
#ax2.yaxis.tick_right()
#ax2.yaxis.set_label_position('right')
ax2.set_ylabel('IPI (s)')
ax2.set_xticklabels([])
ax2.scatter(df_A['year'],df_A['ipi'],s=df_A['peakcounts']*.4,
            linewidths=1.2,
            edgecolors=df_A['stacolors'],c='None')
ax2.scatter(df_B['year'],df_B['ipi'],s=df_B['peakcounts']*.4,
            linewidths=1.2,marker='*',
            edgecolors=df_B['stacolors'],c='None')
ax2.plot(np.array([2002,2014]),np.array([2002,2014])*mod_ipi.params[1] + mod_ipi.params[0],c='black')
ax2.set_xlim([2002,2014])

# Stats
print(mod_ipi.params)
print(mod_ipi.pvalues)


# Subplot 3: Frequency vs year
mod_freq = smf.ols('freq~year',data=df_main).fit()
ax3 = plt.subplot2grid((4,3),(2,0),rowspan=2,colspan=2)
#ax3.yaxis.tick_right()
#ax3.yaxis.set_label_position('right')
ax3.set_ylabel('Frequency (Hz)')
ax3.set_yticklabels(np.append(np.arange(17,25,1),''))
ax3.set_xticks(np.arange(2004,2014,2))
ax3.set_xlim([2002,2014])
ax3.set_xlabel('Year')
ax3.scatter(df_A['year'],df_A['freq'],s=df_A['peakcounts']*.4,
            alpha=1,linewidths=1.2,
            edgecolors=df_A['stacolors'],c='None')
ax3.scatter(df_B['year'],df_B['freq'],s=df_B['peakcounts']*.4,
            alpha=1,linewidths=1.2,marker='*',
            edgecolors=df_B['stacolors'],c='None')
ax3.plot(np.array([2002,2014]),np.array([2002,2014])*mod_freq.params[1] + mod_freq.params[0],c='black')
ax3.set_xlim([2002,2014])

# Stats
print(mod_freq.params)
print(mod_freq.pvalues)


# dummy legend entries
axcol = plt.scatter([],[],s=100,alpha=1,linewidths=1.2,edgecolors=coldict['AX'],c='None',label='Axial')
kenecol = plt.scatter([],[],s=100,alpha=1,linewidths=1.2,edgecolors=coldict['KENE'],c='None',label='KENE')
kemfcol = plt.scatter([],[],s=100,alpha=1,linewidths=1.2,edgecolors=coldict['KEMF'],c='None',label='KEMF')

pcounts = df_decade['propcounts']
pkcounts = df_decade['peakcounts']
# legend circles
s1 = plt.scatter([],[],s=100*.4,edgecolors='grey',c='None',alpha=1,label='100',linewidths=1.2)
s2 = plt.scatter([],[],s=300*.4,edgecolors='grey',c='None',alpha=1,label='300',linewidths=1.2)
s3 = plt.scatter([],[],s=800*.4,edgecolors='grey',c='None',alpha=1,label='800',linewidths=1.2)
# legend stars
s4 = plt.scatter([],[],s=100*.4,edgecolors='grey',c='None',alpha=1,label='100',linewidths=1.2,marker='*')
s5 = plt.scatter([],[],s=300*.4,edgecolors='grey',c='None',alpha=1,label='300',linewidths=1.2,marker='*')
s6 = plt.scatter([],[],s=800*.4,edgecolors='grey',c='None',alpha=1,label='800',linewidths=1.2,marker='*')

legend1 = plt.legend(handles=[axcol,kenecol,kemfcol],loc=3,bbox_to_anchor=(1.09,1.3),scatterpoints=1,frameon=False,fontsize=11,title='Station')
legend2 = plt.legend(handles=[s1,s2,s3],loc=3,bbox_to_anchor=(1.02,0.8),scatterpoints=1,frameon=False,fontsize=11,title='# of seqs: Call A')
legend3 = plt.legend(handles=[s4,s5,s6],loc=3,bbox_to_anchor=(1.02,0.3),scatterpoints=1,frameon=False,fontsize=11,title='# of seqs: Call B')
plt.gca().add_artist(legend1)
plt.gca().add_artist(legend2)


plt.savefig('final-figs/FIG_decade2.tif')



            