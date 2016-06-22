# -*- coding: utf-8 -*-
"""



Created on Sat Jun 18 15:20:29 2016

@author: michw
"""

import numpy as np
import scipy as sp
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
import pandas as pd 
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import copy
import scipy.linalg as LA

#%matplotlib inline# Run this to force plots to open in another window


# Define file list (file name, station, amplitude threshold, inst. name, time span)
flist = [('AX_2006_2007','AX',5,'Ax','06-07'),
         ('AX_2007_2008','AX',5,'Ax','07-08'),
         ('AX_2008_2009','AX',5,'Ax','08-09'),
         ('AX_2009_2010','AX',5,'Ax','09-10'),
         ('AX_2010_2011','AX',5,'Ax','10-11'),
         ('AX_2011_2012','AX',5,'Ax','11-12'),
         ('AX_2012_2013','AX',5,'Ax','12-13')]
         


#flist = [('J63A_2011_2012','CI',10,'J63A'),
#         ('AX_2011_2012','AX',10,'Axial'),
#         ('J23A_2011_2012','CI',10,'J23A'),
#         ('J06A_2011_2012','CI',10,'J06A'),
#         ('G30A_2011_2012','CI',10,'G30A'),
#         ('G03A_2011_2012','CI',10,'G03A')]
#         
         
f1, axarr = plt.subplots(7,5, figsize=(9,9))
axbig = f1.add_subplot(111,frameon=False)
axbig.axes.get_yaxis().set_ticks([])
axbig.axes.get_xaxis().set_ticks([])
#axbig.axes.get_yaxis().set_visible(False)
#axbig.axes.get_xaxis().set_visible(False)

for mdex in np.arange(5):
    plt.setp([a.get_xticklabels() for a in axarr[:-1, mdex]], visible=False)
    
for fdex in np.arange(len(flist)):
    plt.setp([a.get_yticklabels() for a in axarr[fdex, 1:]], visible=False)

#f.tight_layout()
f1.subplots_adjust(hspace=0.05,wspace=0.05)

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
        # Pixel coordinates 
        ipi_pixel_coords = np.arange(len(ipi_centers))       
        freq_pixel_coords = np.arange(len(freq_centers))
        # 1D interpolation to extract pixel coordinates
        f1_ipi = sp.interpolate.interp1d(ipi_pixel_coords,ipi_centers)
        f1_freq = sp.interpolate.interp1d(freq_pixel_coords,freq_centers)
        # Calculate frequency and ipi coordinates for interpolation
        freq_coords = f1_freq(x)
        ipi_coords = f1_ipi(y)
        # Create 2d interp function
        f2 = sp.interpolate.interp2d(freq_centers,ipi_centers,propcount)
        peak_counts = f2(freq_coords,ipi_coords)        
        
        
        plt.clf()
        plt.imshow(np.flipud(f_histo))        
        plt.plot(x,y,'ro')
        
        
        # For straight counts
#        propcount = np.flipud(f_histo)
#        maxval = 70
#        filename = 'monthly2dhists/SMult_AXIAL2.png'
        
        axarr[f,m].scatter(clusters[:,0],clusters[:,1],color='red')                        
        im = axarr[f,m].imshow(propcount,extent=[np.min(freqvec),np.max(freqvec),
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
        axarr[f,m].grid(which='major', axis='both', linestyle='-',color='white')
        
        if m == len(monthlist)-1:
            axarr[f,m].set_ylabel(flist[f][3] + ' ' + flist[f][4], rotation=270, labelpad=13)
            axarr[f,m].yaxis.set_label_position('right')
            
        if f == 0:
            axarr[f,m].set_xlabel(str(monthlist[m]))
            axarr[f,m].xaxis.set_label_position('top')


# Add colorbar
cax,kw = mpl.colorbar.make_axes([ax for ax in axarr.flat])
plt.colorbar(im,cax=cax,extend='max')

# Add big IPI and frequency labels
axbig.set_xlabel('Frequency (Hz)',labelpad=20)
axbig.set_ylabel('IPI (s)',labelpad=25)

#f1.savefig(filename)
