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
import seaborn as sns
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
        allpeaksdat = np.transpose(
            np.array([stavec,datevec,freq_coords,ipi_coords,peak_counts,prop_counts]))
        cols = ['station','datevec','freq','ipi','peakcounts','propcounts']
        
        if ((f==0) & (m==0)): # first instrument/month
            peaksdf = pd.DataFrame(allpeaksdat) # construct data frame
            peaksdf.columns = cols
            #print('thing1')
        else:
            peaksdf_add = pd.DataFrame(allpeaksdat)
            peaksdf_add.columns = cols
            peaksdf = pd.concat([peaksdf,peaksdf_add],axis=0,ignore_index=True) # append to existing data frame
            #print('thing2')
        
#        axarr[f,m].scatter(clusters[:,0],clusters[:,1],color='red') 
        scaling = np.max(peak_counts)/100
        axarr[f,m].scatter(freq_coords,ipi_coords,
                            #s = peak_counts*scaling,
                            #facecolors = 'none', 
                            #edgecolors='red',
                            c='red')   
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

# Convert data types back to numeric 
peaksdf['ipi'] = pd.to_numeric(peaksdf['ipi'])
peaksdf['freq'] = pd.to_numeric(peaksdf['freq'])
peaksdf['peakcounts']= pd.to_numeric(peaksdf['peakcounts'])
peaksdf['propcounts']= pd.to_numeric(peaksdf['propcounts'])
peaksdf['datevec']= pd.to_datetime(peaksdf['datevec'])


plt.figure(5)
plt.clf()
plt.scatter(peaksdf.freq,
            peaksdf.ipi,
            s=peaksdf.propcounts*200,
            alpha=.15,
            edgecolors='none')
plt.grid()

plt.figure(6)
plt.clf()
plt.plot(peaksdf.datevec,
         peaksdf.ipi,'.')
plt.grid()


#f1.savefig('monthly2dhists/SMult_Axial_withdots.png')
#df['date_int'] = df.date.astype(np.int64)