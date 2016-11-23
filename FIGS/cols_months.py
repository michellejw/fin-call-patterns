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
from palettable.colorbrewer.sequential import YlGnBu_6
import pickle

df = pd.read_csv('../SEQ_CODE/ALL_seq_summary.csv')
sta = pd.read_csv('../SEQ_CODE/stationlist.csv')


# Set up colors
unqmonth = [11,12,1,2,3]
colvec = YlGnBu_6.hex_colors[1:]
coldict = {}


for idx in range(5):
    thismo = unqmonth[idx]
    coldict[thismo] = colvec[idx]

# Pickle the dictionary
pickle.dump( coldict, open( "cols_months.p", "wb" ) )

# to load this dictionary in another script (in the same folder): 
#coldict = pickle.load( open( "cols.p", "rb" ) )