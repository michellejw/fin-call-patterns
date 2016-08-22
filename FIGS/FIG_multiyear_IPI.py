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

df = pd.read_csv('../SEQ_CODE/ALL_seq_summary.csv')
sta = pd.read_csv('../SEQ_CODE/stationlist.csv')

df.datevec = pd.to_datetime(df.datevec)
df['year'] = df.datevec.dt.year
df['month'] = df.datevec.dt.month
df['dt'] = df.datevec.dt.date
unqsta = np.unique(df['station'])

# Set up colors (see cols.py)
coldict = pickle.load( open( "cols.p", "rb" ) )    

# loop through each row in df
stacols = []
for idx in range(len(df['station'])):
    stacols.append(coldict[df['station'][idx]])
    
df['stacolors'] = stacols

# Subset the data
df = df[((df.station == 'KENE') | (df.station == 'KEMF') | 
            (df.station=='AX')) & (df.peakcounts > 50)]


#----- PLOTTING -----*
plt.figure(3,figsize=(9,4))
plt.scatter(df['dt'].tolist(),df['ipi'],s=df['peakcounts']*.4,
            edgecolors='black', linewidths=2, c='None',alpha = 0.5)
plt.xlabel('Year',fontsize=12)
plt.ylabel('IPI (s)',fontsize=12)


plt.savefig('FIG_multiyearIPI.tif')




