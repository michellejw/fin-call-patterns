# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 14:54:14 2016

@author: michw
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import pickle
import seaborn as sns
sns.set(style="whitegrid", color_codes=True)
label_size = 12
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

df['songtype'] = df.ipi>=21 # singlet
df['notetype'] = df.freq>=22  # note type B

# Subset the data by date and peakcounts
dfALL = df[(df.peakcounts >50)]
#dfALL = df
df = df[(df.peakcounts > 50) & (df.datevec >= '2011-11-01') & (df.datevec <= '2012-03-01')]
#df = df[(df.datevec >= '2011-11-01') & (df.datevec <= '2012-03-01')]
df = df.sort_values(by='latitude')
# create month-year variable
df['monthyear'] = df['datevec'].apply(lambda x: x.strftime('%b-%Y'))     


dfAsinglet = df[(df.notetype==False) & (df.songtype==True)]
dfAdoublet = df[(df.notetype==False) & (df.songtype==False)]
dfBdoublet = df[(df.notetype==True) & (df.songtype==False)]
dfBsinglet = df[(df.notetype==True) & (df.songtype==True)]
dfA = df[df.notetype==False]
dfB = df[df.notetype==True]


plt.figure(4)
plt.clf()
sns.stripplot(x="monthyear", y="ipi", data=dfAsinglet,jitter=True,marker='o',size=10,alpha=0.5,color="blue");
sns.stripplot(x="monthyear", y="ipi", data=dfAdoublet,jitter=True,marker='o',size=10,alpha=0.5,color="orange");
sns.stripplot(x="monthyear", y="ipi", data=dfBdoublet,jitter=True,marker='o',size=10,alpha=0.5,color="purple");
sns.stripplot(x="monthyear", y="ipi", data=dfBsinglet,jitter=True,marker='o',size=10,alpha=0.5,color="green");
plt.xlabel('Month-Year',fontsize=12)
plt.ylabel('IPI (seconds)',fontsize=12)

# dummy legend entries
Asinglet = plt.scatter([], [],marker='o',s=100,alpha=0.5,color="blue",label='A-singlet');
Adoublet = plt.scatter([], [],marker='o',s=100,alpha=0.5,color="orange",label='A-doublet');
Bdoublet = plt.scatter([], [],marker='o',s=100,alpha=0.5,color="purple",label='B-doublet');
Bsinglet = plt.scatter([], [],marker='o',s=100,alpha=0.5,color="green",label='B-singlet');

plt.legend(handles=[Asinglet,Adoublet,Bsinglet,Bdoublet],loc=2,fontsize=12)
#legend1 = plt.legend(handles=[Asinglet,Adoublet,Bsinglet,Bdoublet],loc=3,bbox_to_anchor=(1.03,1.2),
#                     scatterpoints=1,frameon=False,fontsize=11,title='Note/song type')

plt.savefig('final-figs/FIG_GeoIPI.png')
