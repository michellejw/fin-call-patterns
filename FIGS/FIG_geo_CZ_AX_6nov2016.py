# -*- coding: utf-8 -*-
"""
Created on Sun Nov  6 13:30:58 2016

@author: michw
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import pickle

mpl.rcParams['axes.facecolor'] = 'White'
mpl.rcParams['axes.edgecolor'] = 'gray'
mpl.rcParams['grid.color'] = 'gray'
mpl.rcParams.update({'font.size': 15})


dat = pd.read_csv('CZ_AX_comparison.csv')
sta = pd.read_csv('../SEQ_CODE/stationlist.csv')

# Set up colors (see cols.py)
coldict = pickle.load( open( "cols.p", "rb" ) ) 
# loop through each row in df
stacols = []
for idx in dat.index:
    stacols.append(coldict[dat['station'][idx]])
    
dat['stacolors'] = stacols


dat_y1_CZ = dat[(dat.year==1) & (dat.station=='CZ')]
dat_y1_AX = dat[(dat.year==1) & (dat.station=='AX')]

dat_y2_CZ_SA = dat[(dat.year==2) & (dat.station=='CZ') & (dat.calltype=='SA')]
dat_y2_AX_SA = dat[(dat.year==2) & (dat.station=='AX') & (dat.calltype=='SA')]

dat_y2_CZ_DA = dat[(dat.year==2) & (dat.station=='CZ') & (dat.calltype=='DA')]
dat_y2_AX_DA = dat[(dat.year==2) & (dat.station=='AX') & (dat.calltype=='DA')]

dat_y2_CZ_DB = dat[(dat.year==2) & (dat.station=='CZ') & (dat.calltype=='DB')]
dat_y2_AX_DB = dat[(dat.year==2) & (dat.station=='AX') & (dat.calltype=='DB')]



#----- PLOTTING -----*
plt.figure(1,figsize=(8,8))
plt.clf()
# FREQUENCY
ax1 = plt.subplot2grid((5,3), (0,0),rowspan=2,colspan=2)
# SA, year 1
ax1.scatter(-.05,dat_y1_CZ.freq_mean,c=dat_y1_CZ['stacolors'],s=50,edgecolors='gray')
ax1.scatter(0.05,dat_y1_AX.freq_mean,c=dat_y1_AX['stacolors'],s=50,edgecolors='gray')
# SA, year 2
ax1.scatter(0.95,dat_y2_CZ_SA.freq_mean,c=dat_y2_CZ_SA['stacolors'],s=50,edgecolors='gray',marker='o')
ax1.scatter(1.05,dat_y2_AX_SA.freq_mean,c=dat_y2_AX_SA['stacolors'],s=50,edgecolors='gray',marker='o')
# DA, year 2
ax1.scatter(0.95,dat_y2_CZ_DA.freq_mean,c=dat_y2_CZ_DA['stacolors'],s=50,edgecolors='gray',marker='^')
ax1.scatter(1.05,dat_y2_AX_DA.freq_mean,c=dat_y2_AX_DA['stacolors'],s=50,edgecolors='gray',marker='^')
# DB, year 2
ax1.scatter(0.95,dat_y2_CZ_DB.freq_mean,c=dat_y2_CZ_DB['stacolors'],s=50,edgecolors='gray',marker='s')
ax1.scatter(1.05,dat_y2_AX_DB.freq_mean,c=dat_y2_AX_DB['stacolors'],s=50,edgecolors='gray',marker='s')

ax1.set_xlim([-0.5,1.5])
ax1.set_ylabel('Frequency (Hz)')
plt.grid()
ax1.xaxis.grid(False)

labels = ['2007-2008','2008-2009']
plt.xticks([0,1],labels)
ax1.set_xticklabels([])


# IPI
ax2 = plt.subplot2grid((5,3),(2,0),rowspan=2,colspan=2)
# SA, year 1
ax2.scatter(-0.05,dat_y1_CZ.ipi_mean,c=dat_y1_CZ['stacolors'],s=50,edgecolors='gray')
ax2.scatter(0.05,dat_y1_AX.ipi_mean,c=dat_y1_AX['stacolors'],s=50,edgecolors='gray')
# SA, year 2
ax2.scatter(.95,dat_y2_CZ_SA.ipi_mean,c=dat_y2_CZ_SA['stacolors'],s=50,edgecolors='gray',marker='o')
ax2.scatter(1.05,dat_y2_AX_SA.ipi_mean,c=dat_y2_AX_SA['stacolors'],s=50,edgecolors='gray',marker='o')
# DA, year 2
ax2.scatter(0.95,dat_y2_CZ_DA.ipi_mean,c=dat_y2_CZ_DA['stacolors'],s=50,edgecolors='gray',marker='^')
ax2.scatter(1.05,dat_y2_AX_DA.ipi_mean,c=dat_y2_AX_DA['stacolors'],s=50,edgecolors='gray',marker='^')
# DB, year 2
ax2.scatter(0.95,dat_y2_CZ_DB.ipi_mean,c=dat_y2_CZ_DB['stacolors'],s=50,edgecolors='gray',marker='s')
ax2.scatter(1.05,dat_y2_AX_DB.ipi_mean,c=dat_y2_AX_DB['stacolors'],s=50,edgecolors='gray',marker='s')

ax2.set_xlim([-0.5,1.5])
ax2.set_xticklabels([])
ax2.set_ylabel('IPI (s)')
plt.grid()
ax2.xaxis.grid(False)

labels = ['2007-2008','2008-2009']
plt.xticks([0,1],labels)


# dummy legend entries - colors
AXcol = plt.scatter([],[],s=100,c=coldict['AX'],edgecolors='None',label='AX')
CZcol = plt.scatter([],[],s=100,c=coldict['CZ'],edgecolors='None',label='CZ')

# dummy legend entries - shapes
SA_shp = plt.scatter([],[],s=100,c='gray',edgecolors='None',marker='o',label='Singlet/A')
DA_shp = plt.scatter([],[],s=100,c='gray',edgecolors='None',marker='^',label='Doublet/A')
DB_shp = plt.scatter([],[],s=100,c='gray',edgecolors='None',marker='s',label='Doublet/B')

legend1 = plt.legend(handles=[AXcol,CZcol],loc=3,bbox_to_anchor=(1,1),scatterpoints=1,frameon=False,title='Station',fontsize=13)
legend2 = plt.legend(handles=[SA_shp,DA_shp,DB_shp],loc=3,bbox_to_anchor=(1,0.25),scatterpoints=1,frameon=False,title='Seq./Note type',fontsize=13)
#legend3 = plt.legend(handles=[s4,s5,s6],loc=3,bbox_to_anchor=(1.02,0.3),scatterpoints=1,frameon=False,fontsize=11,title='# of seqs: Call B')
plt.gca().add_artist(legend1)


plt.savefig('final-figs/FIG_AX_CZ_geo.png',dpi=600)