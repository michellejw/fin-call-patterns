# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 19:24:18 2016

@author: michw
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import statsmodels.formula.api as smf
from scipy import stats
import matplotlib as mpl
mpl.rcParams['axes.facecolor'] = 'White'
mpl.rcParams['axes.edgecolor'] = 'lightgray'
mpl.rcParams['grid.color'] = 'lightgray'


dat = pd.read_csv('CZ_AX_comparison.csv')
