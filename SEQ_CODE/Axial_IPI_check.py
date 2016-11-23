# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 08:14:59 2016

@author: michw
"""

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
mpl.rcParams['axes.facecolor'] = 'White'
mpl.rcParams['axes.edgecolor'] = 'darkgray'
mpl.rcParams['grid.color'] = 'darkgray'

#%matplotlib # Run this to force plots to open in another window

# Set up file
flist = [('AX_2012_2013','AX',5,'Axial')]

# Set up plot region       
f1, axarr = plt.subplots(1,1, figsize=(8,8))    

# set vertical offset for line plots
vert = 6
dvert = -0.5

# set histogram bins
binsize = .5
ipi_limits = np.arange(5.,45.,binsize) # bin edges
ipi_centers = ipi_limits[0:-1]+(binsize/2)

monthyear = np.array(((2012,11),(2012,12),
                      (2013,1),(2013,2),(2013,3)),dtype=int)
monthlist = np.array(('Nov','Dec','Jan','Feb','Mar'))



      
f=0
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
    
    # Get histogram values and axes for plotting
    ipihisto = np.histogram(s_ipi,bins=ipi_limits)
    ipihisto = ipihisto[0].astype('float')
    ipihisto2 = np.divide(ipihisto,np.max(ipihisto)) 

    # plot histogram
    axarr.plot(ipi_centers,ipihisto2+vert)
    vert = vert + dvert


axarr.grid()



#
#       
## Set up plot region       
#f1, axarr = plt.subplots(1,1, figsize=(8,8))
#plt.setp([a.get_yticklabels() for a in axarr[:]], visible=False)
#f1.subplots_adjust(hspace=0.05,wspace=0.05)
## remove y axis grid lines
#for a in axarr[:]: plt.setp(a.yaxis.grid(False))
#
## set vertical offset for line plots
#vert = 6
#dvert = -0.5
#
## set histogram bins
#binsize = 1
#ipi_limits = np.arange(5.,45.,binsize) # bin edges
#ipi_centers = ipi_limits[0:-1]+(binsize/2)
#
#monthyear = np.array(((2011,11),(2011,12),
#                      (2012,1),(2012,2),(2012,3)),dtype=int)
#monthlist = np.array(('Nov','Dec','Jan','Feb','Mar'))
#
#         
#for f in np.arange(len(flist)):
#    # loop through each file
#    df = pd.read_hdf(flist[f][0] + '.h5') # read file
#    df['ipi'] = df.ipi  / np.timedelta64(1,'s') # IPI to seconds
#    df = df[(df.isseq==True) & (df.ipi < 60)] # Only keep calls that are in sequence
#    df['m'] = [df.dettime[i].month for i in df.index] # add df column for month
#    df['y'] = [df.dettime[i].year for i in df.index] # add df column for year
#    
#    # Get sequence lengths, assign current seq length to each row of new column seqlen
#    unq_seqs = np.unique(df.seqnum)
#    df['seqlen'] = np.zeros((len(df),1))
#    for s in np.arange(len(unq_seqs)):
#        # get length of current sequence
#        thislen = len(df.loc[(df.seqnum == unq_seqs[s]),'seqnum'])
#        df.loc[(df.seqnum==unq_seqs[s]),('seqlen')] = thislen 
#        
#    df = df[df.seqlen>10]
#    
#    
#    for m in np.arange(len(monthyear)):
#        # subset dataframe containing just the current month in the current file
#        thisgroup = df[(df['m']==monthyear[m][1]) & (df['y']==monthyear[m][0])]
#        sigthresh = np.float(flist[f][2]) # extract threshold defined above for current file
#        hiamps = thisgroup.snr > sigthresh # flag only those calls that exceed amplitude threshold
#        thisgroup = thisgroup[hiamps] # extract calls exceeding ampl. threshold
#        
#        # extract ipi 
#        s_ipi = thisgroup.ipi
#        # find and exclude nans
#        nans = (np.isnan(s_ipi))
#        s_ipi = s_ipi[~nans]
#        
#        # Get histogram values and axes for plotting
#        ipihisto = np.histogram(s_ipi,bins=ipi_limits)
#        ipihisto = ipihisto[0].astype('float')
#        ipihisto2 = np.divide(ipihisto,np.max(ipihisto)) 
#        
#        fig = axarr[m].plot(ipi_centers,ipihisto2+vert)
#        if m == 0:
#            axarr[0].text(-10,vert,flist[f][3])
#
#    vert = vert + dvert
#    
#
#axarr[2].set_xlabel('IPI (s)')      
#for a in np.arange(len(axarr)): 
#    axarr[a].set_xticks(np.arange(10, 45, 5))     
#    axarr[a].set_title(monthlist[a])         
#
#plt.savefig('../FIGS/final-figs/FIG_Geo1Dhistos.png')
