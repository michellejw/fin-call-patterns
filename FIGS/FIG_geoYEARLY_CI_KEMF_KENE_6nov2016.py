# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 15:20:50 2016

@author: michw
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import statsmodels.formula.api as smf
from scipy import stats
import matplotlib as mpl
mpl.rcParams['axes.facecolor'] = 'White'
mpl.rcParams['axes.edgecolor'] = 'gray'
mpl.rcParams['grid.color'] = 'gray'


#df = pd.read_csv('../SEQ_CODE/ALL_seq_summary.csv')
df = pd.read_csv('../SEQ_CODE/ALL_seq_YEARLY_4nov2016_5dBthresh_kurtosis.csv')
#df = pd.read_csv('../SEQ_CODE/ALL_seq_summary_6nov2016_5dBthresh_kurtosis.csv')
sta = pd.read_csv('../SEQ_CODE/stationlist.csv')

df.datevec = pd.to_datetime(df.datevec)
df['year'] = df.datevec.dt.year
df['month'] = df.datevec.dt.month
df['dt'] = df.datevec.dt.date

# Force CI station names to match everything else
df.loc[df.station == 'J23A','station'] = 'CI-J23A'
df.loc[df.station == 'J06A','station'] = 'CI-J06A'
df.loc[df.station == 'J63A','station'] = 'CI-J63A'
df.loc[df.station == 'G03A','station'] = 'CI-G03A'
df.loc[df.station == 'G30A','station'] = 'CI-G30A'

df['stanet'] = df['station']
df = pd.merge(df,sta,how='left',on='stanet')
df['station'] = df['station_x']

# Compute standard error for frequency and ipi measurements
df['ipi_SE'] = df['ipisd']/np.sqrt(df['all_nseq'])
df['freq_SE'] = df['freqsd']/np.sqrt(df['all_nseq'])


# Subset the data            
df = df[(df['datevec'] > '2011-07-01') & (df['datevec'] < '2012-06-30')]

# Set up colors (see cols.py)
coldict = pickle.load( open( "cols.p", "rb" ) )    

# loop through each row in df
stacols = []
for idx in df.index:
    stacols.append(coldict[df['station'][idx]])
    
df['stacolors'] = stacols

# subset stations
df_geo = df[((df.station == 'CI-J23A') | (df.station == 'CI-J06A') | 
            (df.station=='CI-J63A') | (df.station == 'CI-G03A') |
            (df.station=='CI-G30A') | (df.station == 'KEMF') |
            (df.station=='AX')) & (df.peakcounts > 50)]


df_SA = df_geo[(df_geo['ipi'] > 25) & (df_geo['freq'] < 20)]
df_DA = df_geo[(df_geo['ipi'] < 25) & (df_geo['freq'] < 20)]
df_SB = df_geo[(df_geo['ipi'] > 25) & (df_geo['freq'] > 20)]
df_DB = df_geo[(df_geo['ipi'] < 25) & (df_geo['freq'] > 20)]

sumdoublets = np.array(df_DA.ipi) + np.array(df_DB.ipi)

group_SA0 = df_SA.groupby('station').mean()
group_DA0 = df_DA.groupby('station').mean()
group_SB0 = df_SB.groupby('station').mean()
group_DB0 = df_DB.groupby('station').mean()

group_SA = df_SA.groupby('station').first()
group_SA['ipi'] = group_SA0['ipi']

group_DA = df_DA.groupby('station').first()
group_DA['ipi'] = group_DA0['ipi']

group_SB = df_SB.groupby('station').first()
group_SB['ipi'] = group_SB0['ipi']

group_DB = df_SA.groupby('station').first()
group_DB['ipi'] = group_DB0['ipi']

group_sumdoublets = group_DA.ipi + group_DB.ipi


# ----- Plotting ----- #
plt.figure(10)
plt.clf()

ax1 = plt.subplot2grid((1,4),(0,0),colspan=3)
mod_SA = smf.ols('ipi~latitude',data=df_SA).fit()
#plt.errorbar(df_SA['latitude'],df_SA['ipi'],yerr = np.asarray(df_SA['ipi_SE']),fmt='o')
ax1.scatter(df_SA['latitude'],df_SA['ipi'],c=df_SA['stacolors'],s=100,edgecolors='None')
#plt.plot(np.array([39,49]),np.array([39,49])*mod_SA.params[1] + mod_SA.params[0],c='black')

mod_DA = smf.ols('ipi~latitude',data=df_DA).fit()
#plt.errorbar(df_DA['latitude'],df_DA['ipi'],yerr = np.asarray(df_DA['ipi_SE']),fmt='o')
ax1.scatter(df_DA['latitude'],df_DA['ipi'],c=df_DA['stacolors'],s=100,marker='^',edgecolors='None')
#plt.plot(np.array([39,49]),np.array([39,49])*mod_DA.params[1] + mod_DA.params[0],c='black')

mod_DB = smf.ols('ipi~latitude',data=df_DB).fit()
#plt.errorbar(df_DB['latitude'],df_DB['ipi'],yerr = np.asarray(df_DB['ipi_SE']),fmt='o')
ax1.scatter(df_DB['latitude'],df_DB['ipi'],c=df_DB['stacolors'],s=100,marker='s',edgecolors='None')
#plt.plot(np.array([39,49]),np.array([39,49])*mod_DB.params[1] + mod_DB.params[0],c='black')

ax1.scatter(df_DB['latitude'],sumdoublets,c='gray',s=150,marker='*',edgecolors='None',alpha=0.5)

slopeSA = "{:.2f}".format(mod_SA.params[1])
interceptSA = "{:.2f}".format(mod_SA.params[0])
#plt.text(39,30,'y = ' + slopeSA + 'x + ' + interceptSA)

slopeDA = "{:.2f}".format(mod_DA.params[1])
interceptDA = "{:.2f}".format(mod_DA.params[0])
#plt.text(39,20,'y = ' + slopeDA + 'x + ' + interceptDA)

slopeDB = "{:.2f}".format(mod_DB.params[1])
interceptDB = "{:.2f}".format(mod_DB.params[0])
#plt.text(39,12.5,'y = ' + slopeDB + 'x + ' + interceptDB)

plt.ylabel('IPI (s)')
plt.xlabel('Latitude (degrees)')
ax1.grid()

# dummy legend entries - colors
J63Acol = plt.scatter([],[],s=100,c=coldict['CI-J63A'],edgecolors='None',label='CI-J63A')
AXcol = plt.scatter([],[],s=100,c=coldict['AX'],edgecolors='None',label='AX')
KEMFcol = plt.scatter([],[],s=100,c=coldict['KEMF'],edgecolors='None',label='KEMF')
J23Acol = plt.scatter([],[],s=100,c=coldict['CI-J23A'],edgecolors='None',label='CI-J23A')
J06Acol = plt.scatter([],[],s=100,c=coldict['CI-J06A'],edgecolors='None',label='CI-J06A')
G30Acol = plt.scatter([],[],s=100,c=coldict['CI-G30A'],edgecolors='None',label='CI-G30A')
G03Acol = plt.scatter([],[],s=100,c=coldict['CI-G03A'],edgecolors='None',label='CI-G03A')

# dummy legend entries - shapes
SA_shp = plt.scatter([],[],s=100,c='gray',edgecolors='None',marker='o',label='Singlet/A')
DA_shp = plt.scatter([],[],s=100,c='gray',edgecolors='None',marker='^',label='Doublet/A')
DB_shp = plt.scatter([],[],s=100,c='gray',edgecolors='None',marker='s',label='Doublet/B')
sum_shp = plt.scatter([],[],s=150,c='gray',edgecolors='None',marker='*',label='Sum of doublets')

legend1 = plt.legend(handles=[J63Acol,AXcol,KEMFcol,J23Acol,J06Acol,G30Acol,G03Acol],loc=3,bbox_to_anchor=(1,0.5),scatterpoints=1,frameon=False,fontsize=11,title='Station')
legend2 = plt.legend(handles=[SA_shp,DA_shp,DB_shp,sum_shp],loc=3,bbox_to_anchor=(1,0.25),scatterpoints=1,frameon=False,fontsize=11,title='Sequence/Note type')
#legend3 = plt.legend(handles=[s4,s5,s6],loc=3,bbox_to_anchor=(1.02,0.3),scatterpoints=1,frameon=False,fontsize=11,title='# of seqs: Call B')
plt.gca().add_artist(legend1)



# T test to compare means of the higher and lower IPI groups
ipi_hi_SA = df_SA[df_SA['ipi']>30]['ipi']
ipi_lo_SA = df_SA[df_SA['ipi']<=30]['ipi']
stats.ttest_ind(ipi_hi_SA,ipi_lo_SA)



plt.savefig('final-figs/FIG_IPIvLatYEARLY_2011-2012_withsum.png',dpi=600)





