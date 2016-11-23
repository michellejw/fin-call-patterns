# -*- coding: utf-8 -*-
"""
Created on Sun Jun 26 17:50:05 2016

@author: shaun
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

#ool_pds.csv is your data file. Change it as appropriate.
data = pd.read_csv('ool_pds.csv', low_memory=False)

#bug fix for display formats to avoid run time errors
pd.set_option('display.float_format', lambda x:'%f'%x)


#This converts the values in PPETHM, stated Race and Ethnicity, to numeric
data["PPETHM"] = pd.to_numeric(data["PPETHM"], errors='coerce')

##### 
# Why did you convert to numeric? It's already numeric. Unless you sent me a csv file that was already converted 
#####


#This converts the values in W1_C1, stated palitical party affiliation, to numeric
data["W1_C1"] = pd.to_numeric(data["W1_C1"], errors='coerce')
data["W1_C2"] = pd.to_numeric(data["W1_C2"], errors='coerce')


#makes a copy of working data
sub2 = data.copy()

#recode missing values to python missing (NaN) Please Note, PPETHM has no missing/refused data
sub2['W1_C2']=sub2['W1_C2'].replace(-1, np.nan)
sub2['W1_C1']=sub2['W1_C1'].replace(-1, np.nan)


#change format from numeric to categorical as appropriate. So inefficient, considering we made it numerical before.
sub2['PPETHM'] = sub2['PPETHM'].astype('category')
#Rename Categories
sub2['PPETHM']=sub2['PPETHM'].cat.rename_categories(["White", "Black", "Other", "Hispanic", "Multiracial"])

#THIS WORKS AS A WAY TO SHOW COMPARATIVE MEANS BY GROUP
#sns.barplot(x="PPETHM", y="W1_C2", data=sub2)

#Univariate Bar Chart of Ethnicity
#sns.factorplot(sub2["PPETHM"])
#sns.distplot(sub2["PPETHM"], kde=False, norm_hist=True)
#sns.distplot(sub2["PPETHM"])
#plt.xlabel('Stated Ethnicity')
#plt.title('Stated Ethnicity in the Outlook on Life Study')



## THIS WORKS creating Republican vairable 
#def REPUBLICAN (row):
#   if row['W1_C1'] == 1 :
#      return 1
#   elif row['W1_C1'] != 1 :
#      return 0
#
#sub2['REPUBLICAN'] = sub2.apply (lambda row: REPUBLICAN (row),axis=1)
#
###THIS WORKS IF I WANT TO SHOW BIVARIATE ANALYSIS. 
#sns.factorplot(x="PPETHM", y="REPUBLICAN", data=sub2, kind="bar")
#plt.xlabel('Ethnicity')
#plt.ylabel('Republicanism')

#bivariate bar graph C->Q
#sns.factorplot(x='PPETHM', y='W1_C2', data=data, kind="bar", ci=None)
#plt.xlabel('XLAB')
#plt.ylabel('YLAB')


#Creating a QC variable
#print('In order to check the output values, I needed to compare the count of valid cases to the counts in a cross-tabulated table')
#print('This is simply a count of valid cases: the number of people who were not cleaned out when we removed No Answers from W1_C1')
#def ETHNICITYCOUNT (row):
#   if row['W1_C1'] > -1:
#      return 1
#
#sub2['ETHNICITYCOUNT'] = sub2.apply (lambda row: ETHNICITYCOUNT (row),axis=1)
#
#ethdist = sub2['ETHNICITYCOUNT'].value_counts(sort=False, dropna=False)
#print(ethdist)

#crosstabs showing interaction between ethnicity and stated party affiliation
#print('This is a set of counts, showing a cross-tabulation of Stated Ethnicity and Political Party in this dataset')
#print("We are looking for a value of 2245, and that's what is observed.")
#print(pd.crosstab(sub2['PPETHM'], sub2['W1_C1']))
#
#
#print("This the actual desired outcome, a table showing the distribution of stated ethnicity in each stated political party.")
#print("Column Values: 1) Republican. 2) Democrat. 3) Independent. 4) Other.")
#print("Row Values: 1) White, Non-Hispanic. 2) Black, Non-Hispanic. 3) Other, Non-Hispanic. 4) Hispanic. 5) 2+ Races, Non-Hispanic")
#print (pd.crosstab(sub2['PPETHM'], sub2['W1_C1']).apply(lambda r: r/r.sum()))





## Michelle stuff
# Group by party and race
G_party_race = sub2.groupby(['W1_C1','PPETHM']).PPETHM.count()
# Group by party
G_party = sub2.groupby('W1_C1').PPETHM.count()
# Compute percentages of each race within parties
G_pct = G_party_race.div(G_party,level='W1_C1')

# Plot bar chart
ax=G_pct[1].plot.bar() # 1 - I think this corresonds to republican?
ax.set_xlabel('Race',labelpad=10,fontsize=14)
ax.set_ylabel('Ratio of Republicans',labelpad=10,fontsize=14)
plt.setp(ax.get_xticklabels(),rotation=0)

plt.savefig('percentplot.png')




