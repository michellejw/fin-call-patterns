# -*- coding: utf-8 -*-
"""
Create color dictionary for call patterns plots

Created on Mon Aug  1 12:21:53 2016

@author: michw
"""

import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt
#import statsmodels.formula.api as smf
#import statsmodels.api as sm
#from palettable.tableau import GreenOrange_12_r
#from palettable.colorbrewer.qualitative import Set1_9
from palettable.tableau import Tableau_10
import pickle

df = pd.read_csv('../SEQ_CODE/ALL_seq_summary.csv')
sta = pd.read_csv('../SEQ_CODE/stationlist.csv')


# Set up colors
unqsta = np.unique(df['station'])
colvec = Tableau_10.hex_colors
coldict = {}

for idx in range(len(unqsta)):
    thissta = unqsta[idx]
    coldict[thissta] = colvec[idx]

# Pickle the dictionary
pickle.dump( coldict, open( "cols.p", "wb" ) )

# to load this dictionary in another script (in the same folder): 
#coldict = pickle.load( open( "cols.p", "rb" ) )