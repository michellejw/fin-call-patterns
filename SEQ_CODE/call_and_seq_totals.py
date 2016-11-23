# -*- coding: utf-8 -*-
"""
Created on Sat Oct 15 10:01:29 2016

@author: michw
"""

import numpy as np
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
import pandas as pd 
import matplotlib.cm as cm
import copy

#%matplotlib inline# Run this to force plots to open in another window


# Define file list (file name, network name, amplitude threshold, inst. name, time span)
flist = [('AX_2006_2007','AX',5,'AX'),
         ('AX_2007_2008','AX',5,'AX'),
         ('AX_2008_2009','AX',5,'AX'),
         ('AX_2009_2010','AX',5,'AX'),
         ('AX_2010_2011','AX',5,'AX'),
         ('AX_2011_2012','AX',5,'AX'),
         ('AX_2012_2013','AX',5,'AX'),
         ('J63A_2011_2012','CI',5,'J63A'),
         ('J23A_2011_2012','CI',5,'J23A'),
         ('J06A_2011_2012','CI',5,'J06A'),
         ('G30A_2011_2012','CI',5,'G30A'),
         ('G03A_2011_2012','CI',5,'G03A'),
         ('KENE_2003_2004','KECK',12,'KENE'),
         ('KENE_2004_2005','KECK',12,'KENE'),
         ('KENE_2005_2006','KECK',12,'KENE'),
         ('CZ_2007_2008','CZ',10,'CZ'),
         ('CZ_2008_2009','CZ',10,'CZ'),
         ('CZ09_2007_2008','CZ',10,'CZ'),
         ('CZ09_2008_2009','CZ',10,'CZ'),
         ('KEMF_2011_2012','ONC',5,'KEMF'),
         ('KEMF_2012_2013','ONC',5,'KEMF')]
         
# set axis extents
ipibinsize = 1.;
freqbinsize = .4;
ipi_limits = np.arange(5.,45.,ipibinsize)
freq_limits = np.arange(15.,26.,freqbinsize)

# get bin coordinate centers for frequency and IPI
ipi_centers = ipi_limits[:-1] + (ipibinsize/2)
freq_centers = freq_limits[:-1] + (freqbinsize/2)

#monthyear = np.array(((2011,11),(2011,12),
#                      (2012,1),(2012,2),(2012,3)),dtype=int)
thismonth = np.array((11,12,1,2,3),dtype=int)
monthlist = np.array(('Nov','Dec','Jan','Feb','Mar'))

nseqs_ALL = []
ncalls_ALL = []
         
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
    
    # loop through each sequence and label calls w/ sequence length (for filtering short sequences)
    for s in np.arange(len(unq_seqs)):
        # get length of current sequence
        thislen = len(df.loc[(df.seqnum == unq_seqs[s]),'seqnum'])
        df.loc[(df.seqnum==unq_seqs[s]),('seqlen')] = thislen 
    
    # Filter sequences by length
    df = df[df.seqlen>10] 
    # Filter by months (only Nov-Mar were used)
    df = df[(df['m'] >= 11) | (df['m'] <= 3)]
    sigthresh = np.float(flist[f][2])
    hiamps = df.snr > sigthresh
    nans = ((np.isnan(df.ipi))|(np.isnan(df.frequency)))
    ncalls_ALL.append(len(df) - sum(nans))
    nseqs_ALL.append(len(np.unique(df.seqnum)))
    

nseqs_ALL = np.asarray(nseqs_ALL)
ncalls_ALL = np.asarray(ncalls_ALL)

print('Total calls: ' + str(sum(ncalls_ALL)))
print('Total sequences: ' + str(sum(nseqs_ALL)))