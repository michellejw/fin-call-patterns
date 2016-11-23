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
#import matplotlib.pyplot as plt
#from scipy import interpolate
from scipy import stats
import seaborn as sns


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
 
for f in np.arange(len(flist)):
    # loop through each file
    df = pd.read_hdf(flist[f][0] + '.h5') # read file
    df['ipi'] = df.ipi  / np.timedelta64(1,'s') # IPI to seconds
    df = df[(df.isseq==True) & (df.ipi < 60)] # Only keep calls that are in sequence
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
    
    localS = []
    localD = []
    for m in np.arange(len(monthyear)):
        # subset dataframe containing just the current month in the current file
        thisgroup = df[(df['m']==monthyear[m][1]) & (df['y']==monthyear[m][0])]
        sigthresh = np.float(flist[f][2]) # extract threshold defined above for current file
        hiamps = thisgroup.snr > sigthresh # flag only those calls that exceed amplitude threshold
        thisgroup = thisgroup[hiamps] # extract calls exceeding ampl. threshold
        
        # extract ipi 
        s_ipi = thisgroup.ipi
        # find and exclude nans
        nans = (np.isnan(s_ipi))
        s_ipi = s_ipi[~nans]
        singletdata = s_ipi[s_ipi >=22]
        doubletdata = s_ipi[s_ipi <22]
        
        # Storing data for 1 way anova test 
        localS.append(singletdata)
        localD.append(doubletdata)
        

    allS.append(localS)
    allD.append(localD)
    

sstats = np.asarray(allS)
dstats = np.asarray(allD)

                          
# pull data into monthly lists                         
Snov = np.array([sstats[0][0],sstats[1][0],sstats[2][0],sstats[3][0],
                          sstats[4][0],sstats[5][0],sstats[6][0]])

Sdec = np.array([sstats[0][1],sstats[1][1],sstats[2][1],sstats[3][1],
                          sstats[4][1],sstats[5][1],sstats[6][1]])
                          
Sjan = np.array([sstats[0][2],sstats[1][2],sstats[2][2],sstats[3][2],
                          sstats[4][2],sstats[5][2],sstats[6][2]])

Sfeb = np.array([sstats[0][3],sstats[1][3],sstats[2][3],sstats[3][3],
                          sstats[4][3],sstats[5][3],sstats[6][3]])

Smar = np.array([sstats[0][4],sstats[1][4],sstats[2][4],sstats[3][4],
                          sstats[4][4],sstats[5][4],sstats[6][4]])

                          
Dnov = np.array([dstats[0][0],dstats[1][0],dstats[2][0],dstats[3][0],
                          dstats[4][0],dstats[5][0],dstats[6][0]])

Ddec = np.array([dstats[0][1],dstats[1][1],dstats[2][1],dstats[3][1],
                          dstats[4][1],dstats[5][1],dstats[6][1]])
                          
Sjan = np.array([dstats[0][2],dstats[1][2],dstats[2][2],dstats[3][2],
                          dstats[4][2],dstats[5][2],dstats[6][2]])

Sfeb = np.array([dstats[0][3],dstats[1][3],dstats[2][3],dstats[3][3],
                          dstats[4][3],dstats[5][3],dstats[6][3]])

Smar = np.array([dstats[0][4],dstats[1][4],dstats[2][4],dstats[3][4],
                          dstats[4][4],dstats[5][4],dstats[6][4]])
                          

# 1-way ANOVA test to look at whether the ipi peaks for singlets and doublets 
# within a month are statistically similar across stations

Snov_anova = stats.f_oneway(sstats[0][0],sstats[1][0],sstats[2][0],sstats[3][0],
                          sstats[4][0],sstats[5][0],sstats[6][0])

Sdec_anova = stats.f_oneway(sstats[0][1],sstats[1][1],sstats[2][1],sstats[3][1],
                          sstats[4][1],sstats[5][1],sstats[6][1])

Sjan_anova = stats.f_oneway(sstats[0][2],sstats[1][2],sstats[2][2],sstats[3][2],
                          sstats[4][2],sstats[5][2],sstats[6][2])

Sfeb_anova = stats.f_oneway(sstats[0][3],sstats[1][3],sstats[2][3],sstats[3][3],
                          sstats[4][3],sstats[5][3],sstats[6][3])

Smar_anova = stats.f_oneway(sstats[0][4],sstats[1][4],sstats[2][4],sstats[3][4],
                          sstats[4][4],sstats[5][4],sstats[6][4])
 

Dnov_anova = stats.f_oneway(dstats[0][0],dstats[1][0],dstats[2][0],dstats[3][0],
                          dstats[4][0],dstats[5][0],dstats[6][0])

Ddec_anova = stats.f_oneway(dstats[0][1],dstats[1][1],dstats[2][1],dstats[3][1],
                          dstats[4][1],dstats[5][1],dstats[6][1])

Djan_anova = stats.f_oneway(dstats[0][2],dstats[1][2],dstats[2][2],dstats[3][2],
                          dstats[4][2],dstats[5][2],dstats[6][2])

Dfeb_anova = stats.f_oneway(dstats[0][3],dstats[1][3],dstats[2][3],dstats[3][3],
                          dstats[4][3],dstats[5][3],dstats[6][3])

Dmar_anova = stats.f_oneway(dstats[0][4],dstats[1][4],dstats[2][4],dstats[3][4],
                          dstats[4][4],dstats[5][4],dstats[6][4])                         
                          
                          
                          
# Print summary statistics
print('Singlets')
print('Nov J63A = ' + str(np.mean(Snov[0])))
print('Nov AX = ' + str(np.mean(Snov[1])))
print('Nov KEMF = ' + str(np.mean(Snov[2])))
print('Nov J23A = ' + str(np.mean(Snov[3])))
print('Nov J06A = ' + str(np.mean(Snov[4])))
print('Nov G30A = ' + str(np.mean(Snov[5])))
print('Nov G03A = ' + str(np.mean(Snov[6])))








