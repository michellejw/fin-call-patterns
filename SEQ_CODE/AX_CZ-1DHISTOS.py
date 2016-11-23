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
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import copy

#%matplotlib # Run this to force plots to open in another window


# Define file list (file name, station, amplitude threshold, inst. name)
flist = [('AX_2007_2008','AX',5,'Axial'),
         ('AX_2008_2009','AX',5,'Axial'),
         ('CZ_2007_2008','CZ',5,'CZ'),
         ('CZ_2008_2009','CZ',5,'CZ')]

          
# Set up plot region       
f1, axarr = plt.subplots(4,2, figsize=(10,10))
for mdex in np.arange(2):
    plt.setp([a.get_xticklabels() for a in axarr[:-1, mdex]], visible=False)
    
for fdex in np.arange(len(flist)):
    plt.setp([a.get_yticklabels() for a in axarr[fdex, 1:]], visible=False)

#f.tight_layout()
f1.subplots_adjust(hspace=0.05,wspace=0.05)
    


# set vertical offset for line plots
vertstart = 6
dvert = -0.5

# set histogram bins
binsize = .5
ipi_limits = np.arange(10.,40.,binsize) # bin edges
ipi_centers = ipi_limits[0:-1]+(binsize/2) # bin centers for plotting etc

fbinsize = .1
f_limits = np.arange(17.,23.,fbinsize) # bin edges
f_centers = f_limits[0:-1]+(fbinsize/2) # bin centers for plotting etc

months = np.array((11,12,1,2,3),dtype=int)
monthyear = np.array(((2011,11),(2011,12),
                      (2012,1),(2012,2),(2012,3)),dtype=int)
monthlist = np.array(('Nov','Dec','Jan','Feb','Mar'))
      
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
    
    vert = vertstart
    
    for m in np.arange(len(months)):
        # subset dataframe containing just the current month in the current file
        thisgroup = df[(df['m']==monthyear[m][1])]
        sigthresh = np.float(flist[f][2]) # extract threshold defined above for current file
        hiamps = thisgroup.snr > sigthresh # flag only those calls that exceed amplitude threshold
        thisgroup = thisgroup[hiamps] # extract calls exceeding ampl. threshold
    
        # extract ipi and frequency
        s_ipi = thisgroup.ipi
        s_freq = thisgroup.frequency
        # find and exclude nans
        nans = ((np.isnan(s_ipi))|(np.isnan(s_freq)))
        s_ipi = s_ipi[~nans]
        s_freq = s_freq[~nans]    
    
        # Get histogram values and axes for plotting
        ipihisto = np.histogram(s_ipi,bins=ipi_limits)
        ipihisto = ipihisto[0].astype('float')
        ipihisto2 = np.divide(ipihisto,np.max(ipihisto))   
        # Get histogram values and axes for plotting
        fhisto = np.histogram(s_freq,bins=f_limits)
        fhisto = fhisto[0].astype('float')
        fhisto2 = np.divide(fhisto,np.max(fhisto))           

        fig1 = axarr[f,0].plot(ipi_centers,ipihisto2+vert)
        fig2 = axarr[f,1].plot(f_centers,fhisto2+vert)

        vert = vert + dvert




#
#plt.savefig('../FIGS/final-figs/FIG_Geo1Dhistos.png')
