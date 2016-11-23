# -*- coding: utf-8 -*-
"""
Small multiples: Geographic variation

Size: 12x6 panels
2011-2012 - monthly
- Axial
- CI (J63, J23, J06, G30, G03)

Created on Fri May 13 14:19:09 2016

@author: michw
"""

import numpy as np
import pandas as pd 
import pickle
import matplotlib.pyplot as plt
#from scipy import interpolate
from scipy import stats
import seaborn as sns
sns.set_style('whitegrid')
sns.set(font_scale=1.3)
#%matplotlib # Run this to force plots to open in another window


# Define file list (file name, station, amplitude threshold, inst. name)
#flist = [('KENE_2003_2004','KENE',12),
#         ('KENE_2004_2005','KENE',12),
#         ('KENE_2005_2006','KENE',16),
#         ('AX_2006_2007','AX',10),
#         ('AX_2007_2008','AX',10),
#         ('AX_2008_2009','AX',10),
#         ('AX_2009_2010','AX',10),
#         ('AX_2010_2011','AX',10),
#         ('AX_2011_2012','AX',10),
#         ('AX_2012_2013','AX',10)]
         


flist = [('AX_2007_2008','AX',5,'Axial'),
         ('AX_2008_2009','AX',5,'Axial'),
         ('CZ09_2007_2008','CZ',5,'CZ'),
         ('CZ09_2008_2009','CZ',5,'CZ')]
         
         
         
       

monthyear = np.array(((2011,11),(2011,12),
                      (2012,1),(2012,2),(2012,3)),dtype=int)
monthlist = np.array(('Nov','Dec','Jan','Feb','Mar'))

# pre-allocate summary array
dfALL = pd.DataFrame([])
 
for f in np.arange(len(flist)):
    # loop through each file
    df = pd.read_hdf(flist[f][0] + '.h5') # read file
    df['ipi'] = df.ipi  / np.timedelta64(1,'s') # IPI to seconds
    df = df[(df.isseq==True) & (df.ipi < 40)] # Only keep calls that are in sequence
    df['m'] = [df.dettime[i].month for i in df.index] # add df column for month
    df['y'] = [df.dettime[i].year for i in df.index] # add df column for year
    
    # Get sequence lengths, assign current seq length to each row of new column seqlen
    unq_seqs = np.unique(df.seqnum)
    df['seqlen'] = np.zeros((len(df),1))
    for s in np.arange(len(unq_seqs)):
        # get length of current sequence
        thislen = len(df.loc[(df.seqnum == unq_seqs[s]),'seqnum'])
        df.loc[(df.seqnum==unq_seqs[s]),('seqlen')] = thislen 
        
    df = df[df.seqlen>10]
    dfALL = dfALL.append(df)


# Create month-year variable
dfALL['monthyear'] = dfALL['dettime'].apply(lambda x: x.strftime('%b-%Y'))
# Signal amplitude threshold
dfALL = dfALL[dfALL['snr']>5]
dfALL_AX_Y1 = dfALL[(dfALL.dettime >= '2007-11-01') & (dfALL.dettime <= '2008-03-31') & (dfALL['station']=='AX')]
dfALL_CZ_Y1 = dfALL[(dfALL.dettime >= '2007-11-01') & (dfALL.dettime <= '2008-03-31') & (dfALL['station']=='sta09')]
dfALL_AX_Y2 = dfALL[(dfALL.dettime >= '2008-11-01') & (dfALL.dettime <= '2009-03-31') & (dfALL['station']=='AX')]
dfALL_CZ_Y2 = dfALL[(dfALL.dettime >= '2008-11-01') & (dfALL.dettime <= '2009-03-31') & (dfALL['station']=='sta09')]


# Dataframe subsets for singlets and doublets
dfS_nov1 = dfALL[(dfALL['ipi']>=15) & (dfALL['monthyear']=='Nov-2007')]
dfS_dec1 = dfALL[(dfALL['ipi']>=15) & (dfALL['monthyear']=='Dec-2007')]
dfS_jan1 = dfALL[(dfALL['ipi']>=15) & (dfALL['monthyear']=='Jan-2008')]
dfS_feb1 = dfALL[(dfALL['ipi']>=15) & (dfALL['monthyear']=='Feb-2008')]
dfS_mar1 = dfALL[(dfALL['ipi']>=15) & (dfALL['monthyear']=='Mar-2008')]

dfD_nov1 = dfALL[(dfALL['ipi']<15) & (dfALL['monthyear']=='Nov-2007')]
dfD_dec1 = dfALL[(dfALL['ipi']<15) & (dfALL['monthyear']=='Dec-2007')]
dfD_jan1 = dfALL[(dfALL['ipi']<15) & (dfALL['monthyear']=='Jan-2008')]
dfD_feb1 = dfALL[(dfALL['ipi']<15) & (dfALL['monthyear']=='Feb-2008')]
dfD_mar1 = dfALL[(dfALL['ipi']<15) & (dfALL['monthyear']=='Mar-2008')]

dfS_nov2 = dfALL[(dfALL['ipi']>=15) & (dfALL['monthyear']=='Nov-2008')]
dfS_dec2 = dfALL[(dfALL['ipi']>=15) & (dfALL['monthyear']=='Dec-2008')]
dfS_jan2 = dfALL[(dfALL['ipi']>=15) & (dfALL['monthyear']=='Jan-2009')]
dfS_feb2 = dfALL[(dfALL['ipi']>=15) & (dfALL['monthyear']=='Feb-2009')]
dfS_mar2 = dfALL[(dfALL['ipi']>=15) & (dfALL['monthyear']=='Mar-2009')]

dfD_nov2 = dfALL[(dfALL['ipi']<15) & (dfALL['monthyear']=='Nov-2008')]
dfD_dec2 = dfALL[(dfALL['ipi']<15) & (dfALL['monthyear']=='Dec-2008')]
dfD_jan2 = dfALL[(dfALL['ipi']<15) & (dfALL['monthyear']=='Jan-2009')]
dfD_feb2 = dfALL[(dfALL['ipi']<15) & (dfALL['monthyear']=='Feb-2009')]
dfD_mar2 = dfALL[(dfALL['ipi']<15) & (dfALL['monthyear']=='Mar-2009')]


f1,axarr = plt.subplots(2,5,figsize=(12,6))
axbig = f1.add_subplot(111,frameon=False)
axbig.axes.get_yaxis().set_ticks([])
axbig.axes.get_xaxis().set_ticks([])

for mdex in np.arange(5):
    plt.setp([a.get_xticklabels() for a in axarr[:-1, mdex]], visible=False)
    
for ydex in np.arange(2):
    plt.setp([a.get_yticklabels() for a in axarr[ydex, 1:]], visible=False)

f1.subplots_adjust(wspace=0.05,hspace=0.05)

fsize = 0

stacols = sns.color_palette('Blues',2)

sns.boxplot(x='station',y='ipi',data=dfS_nov1,fliersize=fsize,palette = stacols,ax=axarr[0,0],whis=[5,95])
sns.boxplot(x='station',y='ipi',data=dfS_dec1,fliersize=fsize,palette = stacols,ax=axarr[0,1],whis=[5,95])
sns.boxplot(x='station',y='ipi',data=dfS_jan1,fliersize=fsize,palette = stacols,ax=axarr[0,2],whis=[5,95])
sns.boxplot(x='station',y='ipi',data=dfS_feb1,fliersize=fsize,palette = stacols,ax=axarr[0,3],whis=[5,95])
sns.boxplot(x='station',y='ipi',data=dfS_mar1,fliersize=fsize,palette = stacols,ax=axarr[0,4],whis=[5,95])

for m in np.arange(1,5):
    axarr[0,m].set_ylabel('')

for m in np.arange(5):
    axarr[0,m].set_xlabel('')

axarr[0,0].set_ylabel('IPI (s)')

plt.setp([a.get_xticklabels() for a in axarr[1, :]], visible=False)

sns.boxplot(x='station',y='ipi',data=dfS_nov2,fliersize=fsize,palette = stacols,ax=axarr[1,0],whis=[5,95])
sns.boxplot(x='station',y='ipi',data=dfS_dec2,fliersize=fsize,palette = stacols,ax=axarr[1,1],whis=[5,95])
sns.boxplot(x='station',y='ipi',data=dfS_jan2,fliersize=fsize,palette = stacols,ax=axarr[1,2],whis=[5,95])
sns.boxplot(x='station',y='ipi',data=dfS_feb2,fliersize=fsize,palette = stacols,ax=axarr[1,3],whis=[5,95])
sns.boxplot(x='station',y='ipi',data=dfS_mar2,fliersize=fsize,palette = stacols,ax=axarr[1,4],whis=[5,95])

for m in np.arange(1,5):
    axarr[1,m].set_ylabel('')

for m in np.arange(5):
    axarr[1,m].set_xlabel(monthlist[m])

axarr[1,0].set_ylabel('IPI (s)')


# dummy legend entries
Axial = plt.scatter([], [],marker='s',s=100,alpha=1,color=stacols[0],label='Axial');
CZ = plt.scatter([], [],marker='s',s=100,alpha=1,color=stacols[1],label='CZ');

plt.legend(handles=[Axial,CZ],bbox_to_anchor=(1.12,1),fontsize=14,scatterpoints=1)

f1.savefig('../FIGS/final-figs/FIG_boxplots_geo_CZAX.png')
