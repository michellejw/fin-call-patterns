# -*- coding: utf-8 -*-
"""
TESTING
Created on Tue Aug  2 08:51:26 2016

In my final results it looks like there is not much 
frequency variability. Maybe this is correct but I'm 

@author: michw
"""

import numpy as np
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
import pandas as pd 
import matplotlib.cm as cm
import copy
import matplotlib.pyplot as plt


#%matplotlib inline# Run this to force plots to open in another window

# Define file list (file name, network name, amplitude threshold, inst. name, time span)
flist = [('CZ_2007_2008','CZ',5,'CZ'),
         ('CZ_2008_2009','CZ',5,'CZ'),
         ('CZ05_2008_2009','CZ',5,'CZ'),
         ('CZ05_2008_2009','CZ',5,'CZ'),
         ('CZ09_2008_2009','CZ',5,'CZ'),
         ('CZ09_2008_2009','CZ',5,'CZ'),
         ('KEMF_2011_2012','ONC',5,'KEMF'),
         ('KEMF_2012_2013','ONC',5,'KEMF')]
         
f1, axarr = plt.subplots(1,1, figsize=(7,7))
        
         
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
    
    for m in np.arange(len(thismonth)):
        # subset dataframe containing just the current month in the current file
#        thisgroup = df[(df['m']==monthyear[m][1]) & (df['y']==monthyear[m][0])]
        thisgroup = df[(df['m']==thismonth[m])]
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
        
        plt.figure(m+1)
        plt.hist(s_freq,bins=freq_limits)
        
        # Get histogram values and axes for plotting
        f_histo, ipivec, freqvec = np.histogram2d(s_ipi,s_freq,
                                   bins=(ipi_limits,freq_limits))
        my_cmap = copy.copy(cm.get_cmap('viridis'))
        #my_cmap.set_bad('darkgray', alpha=.5)
        my_cmap.set_over('yellow')
        
        totalcount = np.sum(f_histo)
        
        # for proportional counts
        propcount = (np.flipud(f_histo)/totalcount)*100
        maxval = np.max(propcount)*.6
        filename = 'monthly2dhists/SMult_AXIAL3.png'
        
        # For straight counts
#        propcount = np.flipud(f_histo)
#        maxval = 70
#        filename = 'monthly2dhists/SMult_AXIAL2.png'
        
        # Find local maxima
        neighborhood_size = 3;
        threshold = .2;
        data_max = filters.maximum_filter(propcount,neighborhood_size)
        maxima = (propcount==data_max)
        data_min = filters.minimum_filter(propcount,neighborhood_size)
        diff = ((data_max - data_min) > threshold)
        maxima[diff == 0] = 0
        
        labeled, num_objects = ndimage.label(maxima)
        slices = ndimage.find_objects(labeled)
        x,y = [],[]
        for dy,dx in slices:
            x_center = int((dx.start + dx.stop - 1)/2)
            x.append(x_center)
            y_center = int((dy.start + dy.stop - 1)/2)
            y.append(y_center)
            
        # Get call counts at local maxima
        f_histo_flip = np.flipud(f_histo)
        peak_counts = f_histo_flip[y,x]
        prop_counts = propcount[y,x]         
        
        im = axarr.imshow(propcount,extent=[np.min(freqvec),np.max(freqvec),
                                   np.min(ipivec),np.max(ipivec)],
                                   interpolation='nearest',
                                   aspect='auto',
                                   vmin = 0, vmax = maxval,
                                   #norm = LogNorm(),
                                   cmap = my_cmap)
                                   
        totalint = "{:.0f}".format(totalcount);                          
        axarr[f,m].text(18,37,str(totalint),color='white')                           
        axarr[f,m].set_xticks(np.arange(16,25,4))    
        axarr[f,m].set_yticks(np.arange(10,41,10))    

                  