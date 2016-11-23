# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 19:43:42 2016

@author: michw
"""
plt.figure(6)
plt.imshow(propcount,extent=[np.min(freqvec),np.max(freqvec),
       np.min(ipivec),np.max(ipivec)],
       interpolation='nearest',
       aspect='auto',
       vmin = 0, vmax = maxval,
       #norm = LogNorm(),
       cmap = my_cmap)
       
       
plt.plot(freq_coords,ipi_coords,'or')