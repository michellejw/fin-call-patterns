# -*- coding: utf-8 -*-
"""
Created on Sun Nov  6 19:35:17 2016

@author: michw
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
import matplotlib as mpl
from scipy import stats

#import statsmodels.api as sm
#from palettable.tableau import GreenOrange_12
import pickle
mpl.rcParams['axes.facecolor'] = 'White'
mpl.rcParams['axes.edgecolor'] = 'gray'
mpl.rcParams['grid.color'] = 'gray'

#df = pd.read_csv('../SEQ_CODE/ALL_seq_summary.csv')
#df = pd.read_csv('../SEQ_CODE/ALL_seq_summary5b.csv')
df = pd.read_csv('../SEQ_CODE/ALL_seq_summary_6nov2016_5dBthresh_kurtosis.csv')
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




# Set up colors (see cols.py)
coldict = pickle.load( open( "cols.p", "rb" ) )    

# loop through each row in df
stacols = []
for idx in df.index:
    stacols.append(coldict[df['station'][idx]])
    
df['stacolors'] = stacols


# Subset the data            
df1 = df[(df['datevec'] > '2007-10-31') & (df['datevec'] < '2008-04-01')]
df2 = df[(df['datevec'] > '2008-10-31') & (df['datevec'] < '2009-04-01')]

# subset stations
df_geo1 = df1[((df1.station == 'CZ') | (df1.station=='AX')) & (df1.peakcounts > 50)]
df_geo2 = df2[((df2.station == 'CZ') | (df2.station=='AX')) & (df2.peakcounts > 50)]


df_SA1 = df_geo1[(df_geo1['ipi'] >= 25) & (df_geo1['freq'] < 21)]
df_DA1 = df_geo1[(df_geo1['ipi'] < 25) & (df_geo1['freq'] < 21)]
df_SB1 = df_geo1[(df_geo1['ipi'] >= 25) & (df_geo1['freq'] >= 21)]
df_DB1 = df_geo1[(df_geo1['ipi'] < 25) & (df_geo1['freq'] >= 21)]

df_SA2 = df_geo2[(df_geo2['ipi'] >= 25) & (df_geo2['freq'] < 21)]
df_DA2 = df_geo2[(df_geo2['ipi'] < 25) & (df_geo2['freq'] < 21)]
df_SB2 = df_geo2[(df_geo2['ipi'] >= 25) & (df_geo2['freq'] >= 21)]
df_DB2 = df_geo2[(df_geo2['ipi'] < 25) & (df_geo2['freq'] >= 21)]


dat_y1_CZ = df_SA1[(df_SA1.station=='CZ')]
dat_y1_AX = df_SA1[(df_SA1.station=='AX')]

dat_y2_CZ_SA = df_SA2[df_SA2.station=='CZ']
dat_y2_AX_SA = df_SA2[df_SA2.station=='AX']

dat_y2_CZ_DA = df_DA2[df_DA2.station=='CZ']
dat_y2_AX_DA = df_DA2[df_DA2.station=='AX']

dat_y2_CZ_DB = df_DB2[df_DB2.station=='CZ']
dat_y2_AX_DB = df_DB2[df_DB2.station=='AX']



##----- PLOTTING -----*
plt.figure(1,figsize=(8,8))
plt.clf()
# FREQUENCY
ax1 = plt.subplot2grid((5,3), (0,0),rowspan=2,colspan=2)
# SA, year 1
ax1.scatter(-.05,np.mean(dat_y1_CZ.freq),c=dat_y1_CZ['stacolors'],s=50,edgecolors='gray')
ax1.scatter(0.05,np.mean(dat_y1_AX.freq),c=dat_y1_AX['stacolors'],s=50,edgecolors='gray')
# SA, year 2
ax1.scatter(0.95,np.mean(dat_y2_CZ_SA.freq),c=dat_y2_CZ_SA['stacolors'],s=50,edgecolors='gray',marker='o')
ax1.scatter(1.05,np.mean(dat_y2_AX_SA.freq),c=dat_y2_AX_SA['stacolors'],s=50,edgecolors='gray',marker='o')
# DA, year 2
ax1.scatter(0.95,np.mean(dat_y2_CZ_DA.freq),c=dat_y2_CZ_DA['stacolors'],s=50,edgecolors='gray',marker='^')
ax1.scatter(1.05,np.mean(dat_y2_AX_DA.freq),c=dat_y2_AX_DA['stacolors'],s=50,edgecolors='gray',marker='^')
# DB, year 2
ax1.scatter(0.95,np.mean(dat_y2_CZ_DB.freq),c=dat_y2_CZ_DB['stacolors'],s=50,edgecolors='gray',marker='s')
ax1.scatter(1.05,np.mean(dat_y2_AX_DB.freq),c=dat_y2_AX_DB['stacolors'],s=50,edgecolors='gray',marker='s')

ax1.set_xlim([-0.5,1.5])
ax1.set_ylabel('Frequency (Hz)')
plt.grid()
ax1.xaxis.grid(False)

labels = ['2007-2008','2008-2009']
plt.xticks([0,1],labels)
ax1.set_xticklabels([])

# IPI
ax2 = plt.subplot2grid((5,3),(2,0),rowspan=2,colspan=2)
# SA, year 1
ax2.scatter(-0.05,np.mean(dat_y1_CZ.ipi),c=dat_y1_CZ['stacolors'],s=50,edgecolors='gray')
ax2.scatter(0.05,np.mean(dat_y1_AX.ipi),c=dat_y1_AX['stacolors'],s=50,edgecolors='gray')
# SA, year 2
ax2.scatter(.95,np.mean(dat_y2_CZ_SA.ipi),c=dat_y2_CZ_SA['stacolors'],s=50,edgecolors='gray',marker='o')
ax2.scatter(1.05,np.mean(dat_y2_AX_SA.ipi),c=dat_y2_AX_SA['stacolors'],s=50,edgecolors='gray',marker='o')
# DA, year 2
ax2.scatter(0.95,np.mean(dat_y2_CZ_DA.ipi),c=dat_y2_CZ_DA['stacolors'],s=50,edgecolors='gray',marker='^')
ax2.scatter(1.05,np.mean(dat_y2_AX_DA.ipi),c=dat_y2_AX_DA['stacolors'],s=50,edgecolors='gray',marker='^')
# DB, year 2
ax2.scatter(0.95,np.mean(dat_y2_CZ_DB.ipi),c=dat_y2_CZ_DB['stacolors'],s=50,edgecolors='gray',marker='s')
ax2.scatter(1.05,np.mean(dat_y2_AX_DB.ipi),c=dat_y2_AX_DB['stacolors'],s=50,edgecolors='gray',marker='s')

ax2.set_xlim([-0.5,1.5])
ax2.set_xticklabels([])
ax2.set_ylabel('IPI (s)')
plt.grid()
ax2.xaxis.grid(False)

labels = ['2007-2008','2008-2009']
plt.xticks([0,1],labels)


# dummy legend entries - colors
AXcol = plt.scatter([],[],s=100,c=coldict['AX'],edgecolors='None',label='AX')
CZcol = plt.scatter([],[],s=100,c=coldict['CZ'],edgecolors='None',label='CZ')

# dummy legend entries - shapes
SA_shp = plt.scatter([],[],s=100,c='gray',edgecolors='None',marker='o',label='Singlet/A')
DA_shp = plt.scatter([],[],s=100,c='gray',edgecolors='None',marker='^',label='Doublet/A')
DB_shp = plt.scatter([],[],s=100,c='gray',edgecolors='None',marker='s',label='Doublet/B')

legend1 = plt.legend(handles=[AXcol,CZcol],loc=3,bbox_to_anchor=(1,1),scatterpoints=1,frameon=False,title='Station',fontsize=13)
legend2 = plt.legend(handles=[SA_shp,DA_shp,DB_shp],loc=3,bbox_to_anchor=(1,0.25),scatterpoints=1,frameon=False,title='Seq./Note type',fontsize=13)
#legend3 = plt.legend(handles=[s4,s5,s6],loc=3,bbox_to_anchor=(1.02,0.3),scatterpoints=1,frameon=False,fontsize=11,title='# of seqs: Call B')
plt.gca().add_artist(legend1)


plt.savefig('final-figs/FIG_AX_CZ_geo.png',dpi=600)

#two-tailed t-test: frequency
stats.ttest_ind(dat_y1_AX.freq,dat_y1_CZ.freq) # Y1 SA

stats.ttest_ind(dat_y2_AX_SA.freq,dat_y2_CZ_SA.freq) # Y2 SA
stats.ttest_ind(dat_y2_AX_DA.freq,dat_y2_CZ_DA.freq) # Y2 DA
stats.ttest_ind(dat_y2_AX_DB.freq,dat_y2_CZ_DB.freq) # Y2 DB


# two-tailed t-test: ipi
stats.ttest_ind(dat_y1_AX.ipi,dat_y1_CZ.ipi)

stats.ttest_ind(dat_y2_AX_SA.ipi,dat_y2_CZ_SA.ipi)