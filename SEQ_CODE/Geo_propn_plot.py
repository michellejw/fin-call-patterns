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
#import matplotlib.pyplot as plt
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
         


flist = [('J63A_2011_2012','CI',5,'J63A'),
         ('AX_2011_2012','AX',5,'Axial'),
         ('KEMF_2011_2012','ONC',5,'KEMF'),
         ('J23A_2011_2012','CI',5,'J23A'),
         ('J06A_2011_2012','CI',5,'J06A'),
         ('G30A_2011_2012','CI',5,'G30A'),
         ('G03A_2011_2012','CI',5,'G03A')]
         
       

monthyear = np.array(((2011,11),(2011,12),
                      (2012,1),(2012,2),(2012,3)),dtype=int)
monthlist = np.array(('Nov','Dec','Jan','Feb','Mar'))

# pre-allocate stats arrays
means_D = np.empty((len(flist),len(monthlist)))
means_D[:]=np.nan

means_S = np.empty((len(flist),len(monthlist)))
means_S[:] = np.nan        

stdevs_D = np.empty((len(flist),len(monthlist)))
stdevs_D[:]=np.nan

stdevs_S = np.empty((len(flist),len(monthlist)))
stdevs_S[:]=np.nan 

allS = []
allD = []

dfALL = pd.DataFrame([])
 
for f in np.arange(len(flist)):
    # loop through each file
    df = pd.read_hdf(flist[f][0] + '.h5') # read file
    df['ipi'] = df.ipi  / np.timedelta64(1,'s') # IPI to seconds
    df = df[(df.isseq==True) & (df.ipi < 45)] # Only keep calls that are in sequence
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
# Dataframe subsets for each instrument
dfS_nov = dfALL[(dfALL['ipi']>=22) & (dfALL['monthyear']=='Nov-2011')]
dfS_dec = dfALL[(dfALL['ipi']>=22) & (dfALL['monthyear']=='Dec-2011')]
dfS_jan = dfALL[(dfALL['ipi']>=22) & (dfALL['monthyear']=='Jan-2012')]
dfS_feb = dfALL[(dfALL['ipi']>=22) & (dfALL['monthyear']=='Feb-2012')]
dfS_mar = dfALL[(dfALL['ipi']>=22) & (dfALL['monthyear']=='Mar-2012')]

dfD_nov = dfALL[(dfALL['ipi']<22) & (dfALL['monthyear']=='Nov-2011')]
dfD_dec = dfALL[(dfALL['ipi']<22) & (dfALL['monthyear']=='Dec-2011')]
dfD_jan = dfALL[(dfALL['ipi']<22) & (dfALL['monthyear']=='Jan-2012')]
dfD_feb = dfALL[(dfALL['ipi']<22) & (dfALL['monthyear']=='Feb-2012')]
dfD_mar = dfALL[(dfALL['ipi']<22) & (dfALL['monthyear']=='Mar-2012')]



# Set up plot region       
f1, axarr = plt.subplots(5,1, figsize=(12,14))
plt.setp([a.get_xticklabels() for a in axarr[:-1]], visible=False)
f1.subplots_adjust(hspace=0.15,wspace=0.05)


#f1,axarr = plt.subplots(1,5,figsize=(15,3))
#axbig = f1.add_subplot(111,frameon=False)
#axbig.axes.get_yaxis().set_ticks([])
#axbig.axes.get_xaxis().set_ticks([])
#
#plt.setp([a.get_xticklabels() for a in axarr], visible=False)
#plt.setp([a.get_yticklabels() for a in axarr[1:]], visible=False)
#
#f1.subplots_adjust(wspace=0.05)
#
#stacols = sns.color_palette('Blues',7)
#
#sns.boxplot(x='station',y='ipi',data=dfS_nov,fliersize=0,palette = stacols,ax=axarr[0])
#sns.boxplot(x='station',y='ipi',data=dfD_nov,fliersize=0,palette = stacols,ax=axarr[0])
#
#sns.boxplot(x='station',y='ipi',data=dfS_dec,fliersize=0,palette = stacols,ax=axarr[1])
#sns.boxplot(x='station',y='ipi',data=dfD_dec,fliersize=0,palette = stacols,ax=axarr[1])
#
#sns.boxplot(x='station',y='ipi',data=dfS_jan,fliersize=0,palette = stacols,ax=axarr[2])
#sns.boxplot(x='station',y='ipi',data=dfD_jan,fliersize=0,palette = stacols,ax=axarr[2])
#
#sns.boxplot(x='station',y='ipi',data=dfS_feb,fliersize=0,palette = stacols,ax=axarr[3])
#sns.boxplot(x='station',y='ipi',data=dfD_feb,fliersize=0,palette = stacols,ax=axarr[3])
#
#sns.boxplot(x='station',y='ipi',data=dfS_mar,fliersize=0,palette = stacols,ax=axarr[4])
#sns.boxplot(x='station',y='ipi',data=dfD_mar,fliersize=0,palette = stacols,ax=axarr[4])
#
#for m in np.arange(1,len(axarr)):
#    axarr[m].set_ylabel('')
#
#for m in np.arange(len(axarr)):
#    axarr[m].set_xlabel(monthlist[m])
#
#axarr[0].set_ylabel('IPI (s)')
#
## dummy legend entries
#J63A = plt.scatter([], [],marker='s',s=100,alpha=1,color=stacols[0],label='J63A');
#Axial = plt.scatter([], [],marker='s',s=100,alpha=1,color=stacols[1],label='Axial');
#KEMF = plt.scatter([], [],marker='s',s=100,alpha=1,color=stacols[2],label='KEMF');
#J23A = plt.scatter([], [],marker='s',s=100,alpha=1,color=stacols[3],label='J23A');
#J06A = plt.scatter([], [],marker='s',s=100,alpha=1,color=stacols[4],label='J06A');
#G30A = plt.scatter([], [],marker='s',s=100,alpha=1,color=stacols[5],label='G30A');
#G03A = plt.scatter([], [],marker='s',s=100,alpha=1,color=stacols[6],label='G03A');
#
#plt.legend(handles=[J63A,Axial,KEMF,J23A,J06A,G30A,G03A],bbox_to_anchor=(1.12,1),fontsize=14)
#
#f1.savefig('../FIGS/final-figs/FIG_boxplots_geo.png')
