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

df = pd.read_csv('../SEQ_CODE/ALL_seq_summary3.csv')
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
df = df[((df.station == 'KENE') | (df.station == 'KEMF') | 
            (df.station=='AX')) & (df.peakcounts > 50)]
#
#df = df[((df.station == 'KENE') | 
#            (df.station=='AX')) & (df.peakcounts > 50)]

#----- PLOTTING -----*
plt.figure(3,figsize=(9,5))
plt.clf()
plt.scatter(df['dt'].tolist(),df['ipi'],s=40,
            c=df['monthcols'], edgecolors='gray',linewidths=0.5,alpha = 1)
plt.xlabel('Year',fontsize=12)
plt.ylabel('IPI (s)',fontsize=12)

plt.grd()

plt.savefig('final-figs/FIG_multiyearIPI.tif')




