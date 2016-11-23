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
mpl.rcParams['axes.edgecolor'] = 'lightgray'
mpl.rcParams['grid.color'] = 'lightgray'


#df = pd.read_csv('../SEQ_CODE/ALL_seq_summary.csv')
#df = pd.read_csv('../SEQ_CODE/ALL_seq_summary5b.csv')
df = pd.read_csv('../SEQ_CODE/ALL_seq_summary_3nov2016.csv')
sta = pd.read_csv('../SEQ_CODE/stationlist.csv')

df.datevec = pd.to_datetime(df.datevec)
df['year'] = df.datevec.dt.year
df['month'] = df.datevec.dt.month
df['dt'] = df.datevec.dt.date
df['stanet'] = df['station']
df = pd.merge(df,sta,how='left',on='stanet')
df['station'] = df['station_x']

# Compute standard error for frequency and ipi measurements
df['ipi_SE'] = df['ipisd']/np.sqrt(df['all_nseq'])
df['freq_SE'] = df['freqsd']/np.sqrt(df['all_nseq'])

## Set up colors (see cols.py)
#coldict = pickle.load( open( "cols.p", "rb" ) )    
#
## loop through each row in df
#stacols = []
#for idx in df.index:
#    stacols.append(coldict[df['station'][idx]])
#    
#df['stacolors'] = stacols

# Subset the data            
df1 = df[(df['datevec'] > '2007-07-01') & (df['datevec'] < '2008-06-30')]
df2 = df[(df['datevec'] > '2008-07-01') & (df['datevec'] < '2009-06-30')]

# subset stations
df_geo1 = df1[((df1.station == 'AX') | (df1.station == 'CZ')) & (df1.peakcounts > 50)]
df_geo2 = df2[((df2.station == 'AX') | (df2.station == 'CZ')) & (df2.peakcounts > 50)]

df_SA1 = df_geo1[(df_geo1['ipi'] > 25) & (df_geo1['freq'] < 20)]
df_SA2 = df_geo2[(df_geo2['ipi'] > 25) & (df_geo2['freq'] < 20)]

df_DA1 = df_geo1[(df_geo1['ipi'] < 25) & (df_geo1['freq'] < 21)]
df_DA2 = df_geo2[(df_geo2['ipi'] < 25) & (df_geo2['freq'] < 21)]

df_DB1 = df_geo1[(df_geo1['ipi'] < 25) & (df_geo1['freq'] > 21)]
df_DB2 = df_geo2[(df_geo2['ipi'] < 25) & (df_geo2['freq'] > 21)]



# Summarizing IPIs for each station and both years

AX1_ipi = df_SA1[df_SA1.station=='AX']['ipi']
AX1_ipimean = np.mean(AX1_ipi)
AX1_ipisd = np.std(AX1_ipi)
AX2_ipi = df_SA2[df_SA2.station=='AX']['ipi']
AX2_ipimean = np.mean(AX2_ipi)
AX2_ipisd = np.std(AX2_ipi)

CZ1_ipi = df_SA1[df_SA1.station=='CZ']['ipi']
CZ1_ipimean = np.mean(CZ1_ipi)
CZ1_ipisd = np.std(CZ1_ipi)
CZ2_ipi = df_SA2[df_SA2.station=='CZ']['ipi']
CZ2_ipimean = np.mean(CZ2_ipi)
CZ2_ipisd = np.std(CZ2_ipi)


stats.ttest_ind(AX1_ipi,CZ1_ipi)
stats.ttest_ind(AX2_ipi,CZ2_ipi)

###

# Split data into station, year, and call type

# Axial singlet A calls
AX_yr1_SA = df_SA1[df_SA1.station=='AX']
AX_yr2_SA = df_SA2[df_SA2.station=='AX']

# Axial doublet A calls
AX_yr1_DA = df_DA1[df_DA1.station=='AX']
AX_yr2_DA = df_DA2[df_DA2.station=='AX']

# Axial doublet B calls
AX_yr1_DB = df_DB1[df_DB1.station=='AX']
AX_yr2_DB = df_DB2[df_DB2.station=='AX']


# CZ singlet A calls
CZ_yr1_SA = df_SA1[df_SA1.station=='CZ']
CZ_yr2_SA = df_SA2[df_SA2.station=='CZ']

# CZ doublet A calls
CZ_yr1_DA = df_DA1[df_DA1.station=='CZ']
CZ_yr2_DA = df_DA2[(df_DA2.station=='CZ') & (df_DA2.peakcounts > 150)]

# CZ doublet B calls
CZ_yr1_DB = df_DB1[df_DB1.station=='CZ']
CZ_yr2_DB = df_DB2[df_DB2.station=='CZ']





