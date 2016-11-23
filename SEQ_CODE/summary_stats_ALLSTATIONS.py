# -*- coding: utf-8 -*-
"""



Created on Sat Jun 18 15:20:29 2016

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
        
#        # debugging plot
#        plt.clf()
#        plt.imshow(f_histo_flip)        
#        plt.plot(x,y,'ro')                    

        # Calculate frequency and ipi coordinates for local maxima peaks
        freq_coords = freq_centers[x]
        ipi_coords = np.flipud(ipi_centers)[y]
        
        # store freq/ipi peak variables in a pd dataframe
        o = np.ones(freq_coords.shape,dtype=int)
        yr = thisgroup.dettime[thisgroup.index[0]].year
        mo = thisgroup.dettime[thisgroup.index[0]].month
        yearvec = o*yr
        monthvec = o*mo
        datevec = np.tile(pd.to_datetime(yr*10000+mo*100+1,format='%Y%m%d'),len(o))
        stavec = np.chararray(freq_coords.shape,itemsize=4)
        stavec[:] = flist[f][3]
        netvec = np.chararray(freq_coords.shape,itemsize=4)
        netvec[:] = flist[f][1]
        allpeaksdat = np.transpose(
            np.array([stavec,netvec,datevec,freq_coords,ipi_coords,peak_counts,prop_counts]))
        cols = ['station','netvec','datevec','freq','ipi','peakcounts','propcounts']
        
        if ((f==0) & (m==0)): # first instrument/month
            peaksdf = pd.DataFrame(allpeaksdat) # construct data frame
            peaksdf.columns = cols
            #print('thing1')
        else:
            peaksdf_add = pd.DataFrame(allpeaksdat)
            peaksdf_add.columns = cols
            peaksdf = pd.concat([peaksdf,peaksdf_add],axis=0,ignore_index=True) # append to existing data frame
            #print('thing2')
        


# Convert data types back to numeric 
peaksdf['ipi'] = pd.to_numeric(peaksdf['ipi'])
peaksdf['freq'] = pd.to_numeric(peaksdf['freq'])
peaksdf['peakcounts']= pd.to_numeric(peaksdf['peakcounts'])
peaksdf['propcounts']= pd.to_numeric(peaksdf['propcounts'])
peaksdf['datevec']= pd.to_datetime(peaksdf['datevec'])



peaksdf.to_csv('ALL_seq_summary4.csv',sep=',')
