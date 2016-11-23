# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 15:37:01 2016

@author: michw
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#import statsmodels.formula.api as smf
#import statsmodels.api as sm
#from palettable.tableau import GreenOrange_12
import pickle
import matplotlib as mpl

#
#mpl.rcParams['axes.facecolor'] = 'white'
mpl.rcParams['axes.edgecolor'] = 'gray'
mpl.rcParams['grid.color'] = 'gray'
mpl.rcParams.update({'font.size': 15})

#df = pd.read_csv('../SEQ_CODE/ALL_seq_summary_3nov2016.csv')
#df = pd.read_csv('../SEQ_CODE/ALL_seq_summary_4nov2016_12dBthresh.csv')
df = pd.read_csv('../SEQ_CODE/ALL_seq_summary_4nov2016_5dBthresh_kurtosis.csv')
sta = pd.read_csv('../SEQ_CODE/stationlist.csv')

df.datevec = pd.to_datetime(df.datevec)
df['year'] = df.datevec.dt.year
df['month'] = df.datevec.dt.month
df['dt'] = df.datevec.dt.date
unqsta = np.unique(df['station'])

# Set up colors (see cols_months.py)
coldict = pickle.load( open( "cols_months.p", "rb" ) )    

# loop through each row in df
monthcols = []
for idx in range(len(df['month'])):
    monthcols.append(coldict[df['month'][idx]])
    
df['monthcols'] = monthcols

# Subset the data
#df = df[((df.station == 'KENE') | (df.station == 'KEMF') | 
#            (df.station=='AX')) & (df.peakcounts > 50)]
##
df = df[((df.station == 'KENE') | 
            (df.station=='AX')) & (df.peakcounts > 50)]

#----- PLOTTING -----*
plt.figure(3,figsize=(15,5))
plt.clf()
plt.scatter(df['dt'].tolist(),df['ipi'],s=40,
            c=df['monthcols'], edgecolors='gray',linewidths=0.5,alpha = 1)
plt.xlabel('Year',fontsize=15)
plt.ylabel('IPI (s)',fontsize=15)

plt.grid()

# dummy legend entries
Nov = plt.scatter([], [],edgecolors='gray',linewidths=0.5,s=40,alpha=1,color=coldict[11],label='Nov');
Dec = plt.scatter([], [],edgecolors='gray',linewidths=0.5,s=40,alpha=1,color=coldict[12],label='Dec');
Jan = plt.scatter([], [],edgecolors='gray',linewidths=0.5,s=40,alpha=1,color=coldict[1],label='Jan');
Feb = plt.scatter([], [],edgecolors='gray',linewidths=0.5,s=40,alpha=1,color=coldict[2],label='Feb');
Mar = plt.scatter([], [],edgecolors='gray',linewidths=0.5,s=40,alpha=1,color=coldict[3],label='Mar');

plt.legend(handles=[Nov,Dec,Jan,Feb,Mar],bbox_to_anchor=(1.11,1.02),fontsize=15,
           scatterpoints=1,frameon=False)


plt.savefig('final-figs/FIG_multiyearIPI_kcheck.png',dpi=600)




