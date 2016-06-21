# -*- coding: utf-8 -*-
"""
Load Station data

Created on Wed Feb 17 14:38:52 2016

@author: michw
"""

import psycopg2
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import matplotlib
matplotlib.style.use('ggplot')   
from tables import *
#%matplotlib qt


# getdata function
def getdata(tstart,tend,thresh,sta):
    
    conn = psycopg2.connect("""dbname='whaledb' user='postgres' host='localhost'
        password='quakesandwhales'""")
    cur = conn.cursor() 
    querystr = (
        'SELECT findetections.dettime, findetections.frequency,{sep}'
        'findetections.snr, findetections.siglevel, findetections.station {sep}'
        'FROM findetections{sep}'
        'WHERE findetections.dettime > \'{tstart}\'{sep}'
        'AND findetections.dettime <= \'{tend}\'{sep}'
        'AND findetections.snr > {thresh}{sep}'
        'AND findetections.station = \'{sta}\'').format(
            tstart=tstart, tend=tend, thresh=thresh, sta=sta,
            sep='\n')
    cur.execute(querystr)
    the_data = cur.fetchall() 
    colnames = [desc[0] for desc in cur.description]
    the_frame = pd.DataFrame(the_data)
    the_frame.columns = colnames
    cur.close
    conn.close
    return the_frame
    ###
    
## minute-spike correction
#def mspike(df,ipi,basedex):
#     minispikebase = df.dettime[basedex]
#     examplediff = df.dettime-minispikebase
#     exseconds = examplediff / np.timedelta64(1,'s')
#     medmingap = np.median(ipi) / np.timedelta64(1,'s')
#     numsincebase = exseconds/medmingap # minute gaps since baseline
#     mgap_round = np.abs(np.round(numsincebase)-(numsincebase)) # how close to a whole number of minute gaps from the base minute gap?
#     isminspike = mgap_round<1 # true if within 1 second of expected minute spike
#     return isminspike
#     ###
        

# Define stations and years to import - with a file name defined at the end
#datlist = [('2003-07-01 00:00:01','2004-06-30 00:00:01','5','KENE','KENE_2003_2004','1338'),
#           ('2004-07-01 00:00:01','2005-06-30 00:00:01','5','KENE','KENE_2004_2005','12'),
#           ('2005-07-01 00:00:01','2006-06-30 00:00:01','5','KENE','KENE_2005_2006','51')]
           
datlist = [('2003-07-01 00:00:01','2004-06-30 00:00:01','5','KENE','KENE_2003_2004'),
           ('2004-07-01 00:00:01','2005-06-30 00:00:01','5','KENE','KENE_2004_2005'),
           ('2005-07-01 00:00:01','2006-06-30 00:00:01','5','KENE','KENE_2005_2006')]
# NOTE FOR KENE: Added baseline indices for minute-spike removal. These were extracted
# via manual inspection of the data           
           
           
# Loop through datlist to get the data for each instrument/year
for item in datlist:
    tstart0 = item[0]
    tend0 = item[1]
    thresh0 = item[2]
    sta0 = item[3]
    # Call getdata function
    df = getdata(tstart0,tend0,thresh0, sta0)
    df = df.sort_values(by='dettime',ascending=True) # sort by detection time     
    # Add a boolean column to flag calls that are in sequences
    df['isseq'] = np.zeros((len(df),1),dtype=bool)
    df['boutnum'] = np.zeros((len(df),1))
    df['seqnum'] = np.zeros((len(df),1))
    # Add a column of ipi's to the dataframe
    ipi = np.insert(np.diff(df['dettime']),0,0)
        
    # Find and remove minute spikes
    #isminspike = mspike(df,ipi,int(item[5]))   
    ipifloat = ipi / np.timedelta64(1,'s') # convert to float seconds
    isminspike = (ipifloat>58) & (ipifloat<62)
    df = df[~isminspike]
    df['ipi'] = np.insert(np.diff(df['dettime']),0,0)

    # Bout Calculations     
    # Find indices of low IPI calls (less than 40 seconds) 
    lowipidex = df['ipi']<np.timedelta64(40,'s')
    df2 = df[lowipidex]
    # Get a second ipi from the subset of calls with low IPI
    ipi2 = np.insert(np.diff(df.loc[lowipidex,('dettime')]),0,0)
    seqflag = ipi2 > np.timedelta64(20,'m') # find ipi's greater than 20 mins    
    seqflag[0] = True
    seqflag[-1] = True    
    seqdex = df2.index[seqflag] # indices of bout starts
    
    snum = 1   
    
    # Loop through each potential sequence, store sequences containing more than
    # 20 calls with IPIs less than 40 seconds.
    for sdex in np.arange(len(seqdex)-1):
        print('Sequence ' + str(sdex) + ' of ' + str(len(seqdex)-2))
        if sdex == 0:
            sstart = seqdex[0] # bout start
        else:
            sstart = seqdex[sdex]-1 # bout start
        send = seqdex[sdex+1]-1 # bout end
        callsub = (df.index >= sstart) & (df.index <= send) 
        calldex = df.index[callsub]
                
        # Loop through each potential sequence where gaps are less than 20 min.
        # I am keeping sequences of calls where there are at least 20 with 
        # IPI's less than 40 seconds
        if sum(callsub) > 20:
            startseqtime = df.dettime[calldex[0]]
            endseqtime = df.dettime[calldex[-1]]
            #thisseq = df.dettime[(df.dettime >= startseqtime) & (df.dettime <= endseqtime)]
            thisipi = df.ipi[(df.dettime >= startseqtime) & (df.dettime <= endseqtime)] / np.timedelta64(1,'s')
            thisfreq = df.frequency[(df.dettime >= startseqtime) & (df.dettime <= endseqtime)]
            df.loc[(df.dettime >= startseqtime) & (df.dettime <= endseqtime),
                   ('isseq')] = True
            df.loc[(df.dettime >= startseqtime) & (df.dettime <= endseqtime),
                   ('seqnum')] = snum 
            print('sequence length = ' + str(len(thisipi)))
            snum = snum + 1
            
            
    store = pd.HDFStore(item[4]+'.h5')
    store['df'] = df
    store.close()
    
                         
                        
            
            
""" 
    #####
    # Extra stuff...
    #test = df.loc[(df.seqnum == 2),('ipi')]
    ipisec = df[:20]['ipi'] / np.timedelta64(1, 's')
    
    minispikebasedex = 1338
    minispikebase = df.dettime[1338]
    
    #examplediff = df.dettime[0:20]-minispikebase
    #examplediff = df.dettime[-10:]-minispikebase
    examplediff = df.dettime-minispikebase
    exseconds = examplediff / np.timedelta64(1,'s')
    
    medmingap = np.median(df.ipi) / np.timedelta64(1,'s')
    
    numsincebase = exseconds/medmingap # minute gaps since baseline
    
    mgap_round = np.abs(np.round(numsincebase)-(numsincebase)) # how close to a whole number of minute gaps from the base minute gap?
    
    ipifloat = df.ipi / np.timedelta64(1,'s')
    plt.hist(ipifloat,bins=100,range=[0,80])
    
"""