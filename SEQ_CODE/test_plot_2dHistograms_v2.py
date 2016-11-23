# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 08:57:52 2016

@author: michw
"""

import numpy as np
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
import pandas as pd 
import matplotlib.cm as cm
import copy
from scipy.stats import kurtosis
import matplotlib.pyplot as plt


flist = [('AX_2006_2007','AX',5,'AX'), # 0
         ('AX_2007_2008','AX',5,'AX'), # 1
         ('AX_2008_2009','AX',5,'AX'), # 2
         ('AX_2009_2010','AX',5,'AX'), # 3
         ('AX_2010_2011','AX',5,'AX'), # 4
         ('AX_2011_2012','AX',5,'AX'), # 5
         ('AX_2012_2013','AX',5,'AX'), # 6
         ('J63A_2011_2012','CI',5,'J63A'), # 7
         ('J23A_2011_2012','CI',5,'J23A'), # 8
         ('J06A_2011_2012','CI',5,'J06A'), # 9
         ('G30A_2011_2012','CI',5,'G30A'), # 10
         ('G03A_2011_2012','CI',5,'G03A'), # 11
         ('KENE_2003_2004','KECK',12,'KENE'), # 12
         ('KENE_2004_2005','KECK',12,'KENE'), # 13
         ('KENE_2005_2006','KECK',12,'KENE'), # 14
         ('CZ_2007_2008','CZ',5,'CZ'), # 15
         ('CZ_2008_2009','CZ',5,'CZ'), # 16
         ('KEMF_2011_2012','ONC',5,'KEMF'), # 17
         ('KEMF_2012_2013','ONC',5,'KEMF')] # 18
         

# set the file and month   
f = 3
m = 4
         
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





# From the f-loop
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



# From the m-loop
thisgroup = df[(df['m']==thismonth[m])]
sigthresh = np.float(flist[f][2]) # extract threshold defined above for current file
hiamps = thisgroup.snr > sigthresh # flag only those calls that exceed amplitude threshold
thisgroup = thisgroup[hiamps] # extract calls exceeding ampl. threshold
nseqs = len(np.unique(thisgroup.seqnum)) # number of independent sequences - for use in computing stats later
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

# Calculate frequency and ipi coordinates for local maxima peaks
freq_coords = freq_centers[x]
ipi_coords = np.flipud(ipi_centers)[y]



# Plot test image
plt.figure(6)
plt.clf()
plt.imshow(propcount,extent=[np.min(freqvec),np.max(freqvec),
       np.min(ipivec),np.max(ipivec)],
       interpolation='nearest',
       aspect='auto',
       vmin = 0, vmax = maxval,
       #norm = LogNorm(),
       cmap = my_cmap)
       
       
plt.plot(freq_coords,ipi_coords,'or')



# Try picking best peaks using kurtosis (?)

k_ipi = []
k_freq = []

for kdex in np.arange(len(x)):
    peak_ipi = propcount[:,x[kdex]]
    peak_freq = propcount[y[kdex],:]
    k_ipi.append(kurtosis(peak_ipi))
    k_freq.append(kurtosis(peak_freq))

# look for negative kurtosis in either frequency or IPI
k_ipi = np.asarray(k_ipi)
k_freq = np.asarray(k_freq)
ipi_coords2 = ipi_coords[(k_ipi>0) & (k_freq >0) & (peak_counts > 50)]
freq_coords2 = freq_coords[(k_ipi>0) & (k_freq >0) & (peak_counts > 50)]

plt.plot(freq_coords2,ipi_coords2,'kp')
