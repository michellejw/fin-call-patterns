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
from scipy import interpolate
from scipy import stats


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
         
       
# Set up plot region       
f1, axarr = plt.subplots(1,5, figsize=(12,8))
plt.setp([a.get_yticklabels() for a in axarr[:]], visible=False)
f1.subplots_adjust(hspace=0.05,wspace=0.05)
# remove y axis grid lines
for a in axarr[:]: plt.setp(a.yaxis.grid(False))

# set vertical offset for line plots
vert = 6
dvert = -0.5

# set histogram bins
binsize = 1.
ipi_limits = np.arange(5.,45.,binsize) # bin edges
ipi_centers = ipi_limits[0:-1]+(binsize/2)
ipi_centers_D = ipi_centers[ipi_centers<22] #doublet IPI bins
ipi_centers_S = ipi_centers[ipi_centers>=22] # singlet IPI bins

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
        
        
        
        
        
        # Get histogram values and axes for plotting
        ipihisto = np.histogram(s_ipi,bins=ipi_limits)
        ipihisto = ipihisto[0].astype('float')
        ipihisto2 = np.divide(ipihisto,np.max(ipihisto))
        
        # Separate out singlet and doublet peaks
        ipisinglets = ipihisto2[ipi_centers>=22]
        ipidoublets = ipihisto2[ipi_centers<22]
        
        # Linear interpolation: so I can get more precise estimates of 50% half width without doing much math, muahaha
        fd_int = interpolate.interp1d(ipi_centers_D,ipidoublets)
        fs_int = interpolate.interp1d(ipi_centers_S,ipisinglets)
        ipi_centers_Di = np.arange(ipi_centers_D[0],ipi_centers_D[-1],.1)
        ipi_centers_Si = np.arange(ipi_centers_S[0],ipi_centers_S[-1],.1)
        ipidoublets_i = fd_int(ipi_centers_Di)
        ipisinglets_i = fs_int(ipi_centers_Si)
        
        # Find the width of the peak at 50% of maximum amplitude
        s_max = np.max(ipisinglets_i)
        d_max = np.max(ipidoublets_i)
        # places where histogram line crosses the 50% level
        sx_50 = np.sign(np.diff(np.sign(ipisinglets_i-(s_max/2))))
        dx_50 = np.sign(np.diff(np.sign(ipidoublets_i-(d_max/2))))
        # If there is at least one 50% crossing on the left and one 50% crossing on the right
        # SINGLETS 
        if (np.sum(sx_50==1) >= 1) & (np.sum(sx_50==-1) >= 1):
            leftdex = np.where(sx_50==1)[0][0]+1
            rightdex = np.where(sx_50==-1)[0][0]+1
        elif(np.sum(sx_50==1) < 1) & (np.sum(sx_50==-1) >= 1):
            leftdex = 0
            rightdex = np.where(sx_50==-1)[0][0]+1
         elif(np.sum(sx_50==1) >= 1) & (np.sum(sx_50==-1) < 1):
            leftdex = np.where(sx_50==1)[0][0]+1
            rightdex = len(ipisinglets_i)-1
        else:
            leftdex = 0
            rightdex = len(ipisinglets_i)-1
                
        s_width = ipi_centers_Si[rightdex] - ipi_centers_Si[leftdex]
        
        # Find the weighted mean IPI for singlets and for doublets
        s_mean = np.sum(ipisinglets_i[leftdex:rightdex]*ipi_centers_Si[leftdex:rightdex])/np.sum(ipisinglets_i[leftdex:rightdex])
        
        # DOUBLETS
        if (np.sum(dx_50==1) >= 1) & (np.sum(dx_50==-1) >= 1):
            leftdex = np.where(dx_50==1)[0][0]+1
            rightdex = np.where(dx_50==-1)[0][0]+1
        elif(np.sum(dx_50==1) < 1) & (np.sum(dx_50==-1) >= 1):
            leftdex = 0
            rightdex = np.where(dx_50==-1)[0][0]+1
        elif(np.sum(dx_50==1) >= 1) & (np.sum(dx_50==-1) < 1):
            leftdex = np.where(dx_50==1)[0][0]+1
            rightdex = len(ipidoublets_i)-1
        else: 
            leftdex = 0
            rightdex = len(ipidoublets_i)-1
                
        d_width = ipi_centers_Di[rightdex] - ipi_centers_Di[leftdex]
        
        # Find the weighted mean IPI for singlets and for doublets
        d_mean = np.sum(ipidoublets_i[leftdex:rightdex]*ipi_centers_Di[leftdex:rightdex])/np.sum(ipidoublets_i[leftdex:rightdex])
    
        # Put mean and half width of dist at 1/2 peak amplitude
        means_D[f,m] = d_mean
        stdevs_D[f,m] = d_width/2
        means_S[f,m] = s_mean
        stdevs_S[f,m] = s_width/2
        
          
        fig = axarr[m].plot(ipi_centers,ipihisto2+vert)
        if m == 0:
            axarr[0].text(-10,vert,flist[f][3])

        
        
    vert = vert + dvert
    allS.append(localS)
    allD.append(localD)
    

axarr[2].set_xlabel('IPI (s)')      
for a in np.arange(len(axarr)): 
    axarr[a].set_xticks(np.arange(10, 45, 5))     
    axarr[a].set_title(monthlist[a])         

plt.savefig('../FIGS/final-figs/FIG_Geo1Dhistos.png')



# 1-way ANOVA test to look at whether the ipi peaks for singlets and doublets 
# within a month are statistically similar across stations





