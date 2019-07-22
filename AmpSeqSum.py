#!/usr/bin/python

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
# Author:	Todd D. Yoder
#
# Date:	May 2019
#
# Purpose:
# 1. Calculates Indel%
# 2. Concatenates Knock Out and Genotype Composition summary files from AmpSeqQC.py
# 3. Slices DF and transposes individual DFs
# 4. Reformats table


# #Modules:--------------------------------------------------------------------------------------------------------------

import sys
import numpy as np
import pandas as pd
import os
import subprocess
from subprocess import call


# # User Variables:--------------------------------------------------------------------------------------------------------------

path = sys.argv[1]
os.chdir(path)


sample_dict= {"1":"VP-19-002.02.G1/TDN,R1C1 S1D5","2":"VP-19-002.02.G1/TDN,R1C1 S1D7","3":"VP-19-002.02.G1/TDN,R1C1 Post-H","4":"VP-19-002.01.G1 D(-2)",
              "5":"VP-19-002.02.G1/Pre-TDN,D0"}


KO_cat = []
indel_dict = {}

# #__________________Pivot Summary Formating____________________________#
script="cat *pivot_summary.csv>merge.csv"
call(script,shell=True)

pivot_sum = pd.read_csv('merge.csv', sep=",", dtype =object)
pivot_sum = pivot_sum[~pivot_sum['Sample_Name'].str.contains('Sample_Name')]
Sample_list = pivot_sum['Sample_Name'].unique()

pivot_sum2 = pivot_sum[pivot_sum['Variant'].str.contains('D|I')].copy()
pivot_sum2['Percentage'] = pivot_sum2['Percentage'].astype(float)
pivot_sum2= pivot_sum2.pivot_table(index=['Sample_Name'], values=['Percentage'],aggfunc={'Percentage':'sum'})
pivot_sum3 = pd.DataFrame(pivot_sum2.to_records())

pivot_sum_noedit = pivot_sum[pivot_sum['Variant'].str.contains('NoEdit')] #<------- Assumes that no sample is 100% editted***
pivot_sum_noedit= pivot_sum_noedit.pivot_table(index=['Sample_Name'], values=['Percentage'],aggfunc={'Percentage':'sum'})
pivot_sum_noedit2 = pd.DataFrame(pivot_sum_noedit.to_records())
pivot_sum_noedit2['Percentage'] = 0

combo_pi = pd.concat([pivot_sum3,pivot_sum_noedit2],axis=0)
combo_pi_pivot= combo_pi.pivot_table(index=['Sample_Name'], values=['Percentage'],aggfunc={'Percentage':'sum'})
combo_pi_pivot = pd.DataFrame(combo_pi_pivot.to_records())
combo_pi_pivot['Percentage'] = round(combo_pi_pivot['Percentage'],2)


combo_pi_pivot['Samp_No'] = combo_pi_pivot['Sample_Name'].str.split('_').str.get(2)
combo_pi_pivot['Locus'] = combo_pi_pivot['Sample_Name'].str.split('_').str.get(1)

indel_dict = dict(zip(combo_pi_pivot['Sample_Name'], combo_pi_pivot['Percentage']))


# #__________________KO Summary Formating____________________________#

script="cat *_summary_knockout.csv>ko_merge.csv"
call(script,shell=True)

ko_sum = pd.read_csv('ko_merge.csv', sep=",", dtype =object,header=None)

# print(len(sample_dict))

hope = np.array_split(ko_sum, (len(sample_dict)*4)) #<---- 120 rows in merged file / 5 lines per samples => 24

for i in hope:

    i_tran= i.transpose()
    i_tran.columns = ['Sample Name:','Total Read Pairs:','Proper Read Pairs:','Total usable percentage:','Sample Edit Level:']
    KO_cat.append(i_tran) ### appends sub-dfs


KO_df = pd.concat(KO_cat,axis=0)     ### concats to one DF
KO_df= KO_df[KO_df['Sample Name:'].str.contains('19')]#<-------- Assumes all Specimen Names contain 19 from '2019'

KO_df['Samp_No'] = KO_df['Sample Name:'].str.split('_').str.get(2)
KO_df['Locus'] = KO_df['Sample Name:'].str.split('_').str.get(1)

KO_df = KO_df.sort_values(['Samp_No','Locus'], ascending=[True,True])

KO_df["Insertion/Deletion_Level_(%)"] = KO_df['Sample Name:'].map(indel_dict)

KO_df["Specimen_Name"] = KO_df['Samp_No'].map(sample_dict)

KO_df = KO_df[['Specimen_Name','Locus','Total Read Pairs:','Proper Read Pairs:','Total usable percentage:','Sample Edit Level:',"Insertion/Deletion_Level_(%)"]]

KO_df['Total usable percentage:'] = round(KO_df['Total usable percentage:'].astype(float),2)
KO_df['Sample Edit Level:'] = round(KO_df['Sample Edit Level:'].astype(float),2)

# #_______________Output Block_____________________________________________# #
KO_df.to_csv('VP-19-002.02.G1_summary_table.csv', sep=',', header=True,index=False)
# #
print(KO_df)
