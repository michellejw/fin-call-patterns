# -*- coding: utf-8 -*-
"""
CreateMap.py
Build a map that shows my study area, and all instruments, 
overlaid on bathymetry.

Created on Wed May 11 11:29:31 2016
Modified on Mon Aug 1 2016

@author: michw
"""

import numpy as np
import pandas as pd 
from scipy import interpolate
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from palettable.colorbrewer.qualitative import Set1_5

#%matplotlib

minlat = 39
maxlat = 50
minlon = -131
maxlon = -120

bathyfile = 'bathy/callpatternsbathy2.csv'
bathydat = pd.read_csv(bathyfile, low_memory=False)

lat = np.array(bathydat['latitude'])
lon = np.array(bathydat['longitude'])
topo = np.array(bathydat['topo'],dtype = 'float64')

resolution = 0.008333333*2 # http://coastwatch.pfeg.noaa.gov/erddap/griddap/usgsCeSrtm30v6.html

# Determine the number of grid points in the x and y directions
nx = complex(0,(max(lon)-min(lon))/resolution)
ny = complex(0,(max(lat)-min(lat))/resolution)

# Build 2 grids: One with lats and the other with lons
grid_x, grid_y = np.mgrid[min(lon):max(lon):nx,min(lat):max(lat):ny]
 
# Interpolate topo into a grid (x by y dimesions)
grid_z = interpolate.griddata((lon,lat),topo,(grid_x,grid_y),method='linear')
 
# Make an empty 'dictionary'... place the 3 grids in it.
TOPO = {}
TOPO['lats']=grid_y
TOPO['lons']=grid_x
TOPO['topo']=grid_z

 
# Create map
m = Basemap(projection='mill', llcrnrlat=minlat,urcrnrlat=maxlat,llcrnrlon=minlon, urcrnrlon=maxlon,resolution='h')
#x,y = m(grid_x[1:],grid_y[1:])
x,y = m(grid_x,grid_y)
 
fig1 = plt.figure(4, figsize = (6,8))
fig1.clf()
ax1 = fig1.add_subplot(111)
cs = m.pcolor(x,y,grid_z,cmap=plt.cm.Greys_r,vmax=0)
m.fillcontinents(color='tan', lake_color='lightgrey')
m.drawcoastlines(color='darkslategray')
# Add parallels and meridians
parallels = np.arange(40,51,2)
meridians = np.arange(-130,-119,2)
m.drawparallels(parallels,labels=[True, False, False, False],color='lightgray')
m.drawmeridians(meridians,labels=[False,False, True, False],color='lightgray')

stations = pd.read_csv('stationlist_offset.csv')
# stations list without KENE, so we don't over-plot
#stations = stations[stations['station']!='KENE']
net = stations['network']
stlat = np.array(stations['latitude'])
stlon = np.array(stations['longitude'])
depth = np.array(stations['depth'])
stname = np.array(stations['station'])

xsta, ysta = m(stlon, stlat) # get positions in map coordinates

#colvec = GreenOrange_12.hex_colors
colvec = Set1_5.hex_colors

pts_CI = m.scatter(xsta[np.array(net=='CI')],ysta[np.array(net=='CI')],
       s=80,color=colvec[0], edgecolor='darkslategray',label='CI')

pts_CI = m.scatter(xsta[np.array(net=='Keck')],ysta[np.array(net=='Keck')],
       s=80,color=colvec[1], edgecolor='darkslategray',label='Keck')    
       
pts_CI = m.scatter(xsta[np.array(net=='ONC')],ysta[np.array(net=='ONC')],
       s=80,color=colvec[2], edgecolor='darkslategray',label='ONC')   
       
pts_CI = m.scatter(xsta[np.array(net=='Axial')],ysta[np.array(net=='Axial')],
       s=80,color=colvec[3], edgecolor='darkslategray',label='Axial')      

pts_CI = m.scatter(xsta[np.array(net=='Colza')],ysta[np.array(net=='Colza')],
       s=80,color=colvec[4], edgecolor='darkslategray',label='Colza')
       
# Add labels to points
for an in np.arange(len(stations)):
    if (stname[an] != 'KEMF') & (stname[an] != 'KENE'):
        lab = ax1.text(xsta[an],ysta[an]+28000,
                     stname[an],color='darkslategray',
                     fontsize=12)        
    else:
        lab = ax1.text(xsta[an]+28000,ysta[an]-10000,
                     stname[an],color='darkslategray',
                     fontsize=12)           
    
       
       
#ax1.legend(scatterpoints=1)

fig1.savefig('FIG_map.tif')
#fig1.savefig('CallPatterns_map2.eps') # takes forever to open using illustrator