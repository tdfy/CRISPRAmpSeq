#!/usr/bin/python

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
# Author:	Todd D. Yoder
#
# Date:	April 2019
#
# Purpose:
# 1. Generate lists of control and treatment files
# 2. Combine all CIGAR columns from control files into one dataframe
# 3. Filter Genotypes that occur <0.1%
# 4. Call intersections between control files and export list of putative PCR/Sequencing artifacts
# 5. Iterate through control/treatment files and remove putative PCR/sequencing artifacts
# 6. Strip header from input files
# 7. Recalculate metrics for new header
# 8. Concatenate new header on to files w/o technical artifacts
# 9. Parse CIGAR strings containing variants into individual variants with coorelating physical position {Figure Generation}
# 10. Generate 'Edit Type' genotype composition per sample {Figure Generation}

# Note from Kyle:
# Consider as an exercise: write down the major steps in the process, like 5-10 lines
# Then imagine writing a function, one for each step, that does the work and returns a result
# can you scoop up part of your code into one of these functions? If so, this is progress.

# #Modules/Libraries:--------------------------------------------------------------------------------------------------------------

import sys
import numpy as np
import pandas as pd
import glob
import os
from cigar import Cigar
import collections
import argparse

# # User Variables:--------------------------------------------------------------------------------------------------------------
#
#
# # Global Declarations:--------------------------------------------------------------------------------------------------------------
#
CIGAR_cat = []
data = {}
header_dict = {}
pivot_dict = {}
summary = pd.DataFrame()
pivot_summary = pd.DataFrame()
String_list = ['I','D','X']
TRAC3_SNPlist = ['198=1X36=','1S198=1X36=','197=1X36=','198=1X35=1X1S','188=1X46=']
CIITA_SNPlist = ['52=1X176=2S','1X51=11X176=2S','1X51=1X176=2S']
PIK_SNPlist = ['39=1D182=']
B2M_SNPlist = []
#
# print("Control Samples:",control_list)
# print("Test Samples:",filelist)
#




path = sys.argv[1]
os.chdir(path)
config = 'config.csv'
configed = pd.read_csv(path +'/'+ config, sep=",", dtype =object)

Targs = configed['Target'].unique()

for comp in Targs:
    comp = configed.loc[(configed['Target'] == comp)]
    comp = comp.sort_values(['Design'], ascending=[False])

    for samp in comp['Directory']+comp['File_Name']:
        print(comp)



print(samp)

# # :--------------------------------------------------------------------------------------------------------------
# Locus = sys.argv[1]
# # path = sys.argv[3]
# # always assume filepaths are relative to current directory
# # os.chdir(path)
# Run_ID = sys.argv[2]
# Noise = sys.argv[4]
# # filelist = glob.glob(Run_ID+r'_'+Locus+r'_*.output.csv')
# control_list = [Run_ID+r'_'+Locus+r'_4']
# Treatment_list = [Run_ID+r'_'+Locus+r'_5',Run_ID+r'_'+Locus+r'_1',Run_ID+r'_'+Locus+r'_2',Run_ID+r'_'+ Locus+r'_3']
# extension = '.output.csv'
# list_1 = Treatment_list + control_list
# filelist = [x + extension for x in list_1]
#
#
# # Global Declarations:--------------------------------------------------------------------------------------------------------------
#
# CIGAR_cat = []
# data = {}
# header_dict = {}
# pivot_dict = {}
# summary = pd.DataFrame()
# pivot_summary = pd.DataFrame()
# String_list = ['I','D','X']
# TRAC3_SNPlist = ['198=1X36=','1S198=1X36=','197=1X36=','198=1X35=1X1S','188=1X46=']
# CIITA_SNPlist = ['52=1X176=2S','1X51=11X176=2S','1X51=1X176=2S']
# PIK_SNPlist = ['39=1D182=']
# B2M_SNPlist = []
#
# print("Control Samples:",control_list)
# print("Test Samples:",filelist)
#
# # # _________Generate DF_________________________________________________________________________________#
#
# for sample in filelist:
#     data[sample[:-11]] = pd.read_csv(sample, sep=",", dtype =object,skiprows=4)
# #
#     # data[sample[:-11]] = data[sample[:-11]][data[sample[:-11]]['Percentage'].astype(float) >= .1] #<--------- threshold logic
#
#     #####____NEW LOGIC BLOCK____20190710____####
#
#     data[sample[:-11]] = data[sample[:-11]][data[sample[:-11]]['Percentage'].astype(float) >= .01] #<--------- threshold logic (.01)
#
#     ###_________________________________________________________####
#
#
#     # This is a pattern that suggests using a dictionary
#     # snp_lists = {"TRAC": TRAC3_SNPlist, "CIITA": etc... }
#     # snps = snp_lists[Locus]
#     # If the Locus isn't found here, Python will give a nice error
#     if Locus == "TRAC":
#        data[sample[:-11]].iloc[:,4] = np.where((data[sample[:-11]].iloc[:,5].isin(TRAC3_SNPlist)), 'NoEdit', data[sample[:-11]].iloc[:,4])
#     elif Locus == "CIITA":
#         data[sample[:-11]].iloc[:,4] = np.where((data[sample[:-11]].iloc[:,5].isin(CIITA_SNPlist)), 'NoEdit', data[sample[:-11]].iloc[:,4])
#     elif Locus == "PIK":
#         data[sample[:-11]].iloc[:,4] = np.where((data[sample[:-11]].iloc[:,5].isin(PIK_SNPlist)), 'NoEdit', data[sample[:-11]].iloc[:,4])
#     elif Locus == "B2M":
#         data[sample[:-11]].iloc[:,4] = np.where((data[sample[:-11]].iloc[:,5].isin(B2M_SNPlist)), 'NoEdit', data[sample[:-11]].iloc[:,4])
#
#     else:
#         print("Target Error")
#
# for x in [data[x] for x in control_list]:
#     x = x.iloc[:,5]
#     CIGAR_frame =x.to_frame()     ### converts series to df in order to concat into one DF
#     CIGAR_cat.append(CIGAR_frame)     ### appends sub-dfs
#
# CIGAR_df = pd.concat(CIGAR_cat,axis=1,)     ### concats to one DF
# # print(CIGAR_df)
#
# # # ___Call Intersections___#
#
# if len(CIGAR_df.columns) == 2:
#     u = set.intersection(set(CIGAR_df.iloc[:,0]),set(CIGAR_df.iloc[:,1]))
#     print("Number of Controls:",len(CIGAR_df.columns))
#
# elif len(CIGAR_df.columns) == 3:
#         u1 = set.intersection(set(CIGAR_df.iloc[:,0]),set(CIGAR_df.iloc[:,1]))
#         u2 = set.intersection(set(CIGAR_df.iloc[:,0]),set(CIGAR_df.iloc[:,2]))
#         u3 = set.intersection(set(CIGAR_df.iloc[:,1]),set(CIGAR_df.iloc[:,2]))
#         print("Number of Controls:",len(CIGAR_df.columns))
#
#         u= u1|u2|u3
#         print("Length of Artifacts:",len(u1),len(u2),len(u3))
#
# elif len(CIGAR_df.columns) == 1:
#             u = []
#             cig_df = pd.DataFrame()
#             pos_list = []
#             var_list = []
#             cig_pack = []
#
#             changes = CIGAR_df.iloc[:,0].tolist()
#             for cig in changes:
#                 CIGAR_edit = Cigar(cig)
#                 CIG_list = list(CIGAR_edit.items())
#
#                 for index, tup in enumerate(CIG_list):
#                     ch = CIG_list[index][1]
#                     if ch == 'X':
#
#                         slice_l = CIG_list[:index + 1]
#                         pos = str(sum([t[0] for t in slice_l]))
#                         variant = slice_l[-1][1]
#                         pos_list.append(pos)
#
#                     else:
#                         pass
#
#             artifact_pos = [item for item, count in collections.Counter(pos_list).items() if count > 50] #<--------- Threshold for substitution pileup (285) ORG--500
#
#             print(artifact_pos)
#
#             new_list = [int(x)-1 for x in artifact_pos]
#
#             for keys in data:
#                 data[keys]['sus_artifact'] = 'real'
#                 for row in keys:
#                     for position in new_list:
#                         data[keys]['Type'] = np.where((data[keys]['Type']=='S') & (position in data[keys]['CIGAR']),'NoEdit',data[keys]['Type'])
#
#
# else:
#     print("No Control Present")
#     exit()
#
# print("Number of Artifacts:",len(u))
# #
# # # _______Export_list________________________#
#
# with open(Locus+r"_artifacts.txt", "w") as output:
#     output.write(str(u))
#
# # # ___Loop through All Samples in Filelist___#
#
# for file in filelist:
#     with open(file) as myfile:
#         header_dict[file[:-11]] = [next(myfile) for x in range(4)]
#     # print(header_dict[file[:-11]]) #<---------------- HEADER CHECK
#
# # # __________________________QC Analysis_____________________________________#
#
# for keys in data:
#     cig_df = pd.DataFrame()
#     pos_list = []
#     var_list = []
#     cig_pack = []
#     varlen_list = []
#
#     QC = pd.read_csv(file, sep=",", dtype =object,skiprows=4)
#
#     # data[keys]['sus_artifact'] = np.where((data[keys].iloc[:,5].isin(u)) & (data[keys].iloc[:,4] != 'NoEdit')  |
#     #     (data[keys].iloc[:,5] == '*') ,'artifact','real')
#
#
#     # data[keys]['sus_artifact'] = np.where((data[keys].iloc[:,5].str.contains('D|I')==False) & (data[keys].iloc[:,2].astype(float) > 4) &
#     #     (data[keys].iloc[:,4] != 'NoEdit'),'artifact',data[keys]['sus_artifact']) #<----------------- Removes substitution artifacts that do not contain Indel and are observed at a >4% freq
#
#         #####____NEW LOGIC BLOCK____20190710____####
#
#     data[keys]['sus_artifact'] = 'real' ###<<<------------- MUST REVIEW AND WRITE IN NEW LOGIC!
#
#     # data[keys].iloc[:,4] = np.where((data[keys].iloc[:,5].isin(u)) & (data[keys].iloc[:,4] != 'NoEdit')  |
#     #     (data[keys].iloc[:,5] == '*') ,'NoEdit',data[keys].iloc[:,4])
#     #
#     #
#     # data[keys].iloc[:,4] = np.where((data[keys].iloc[:,5].str.contains('D|I')==False),'NoEdit', data[keys].iloc[:,4]) #<----------------- Redefines reads that are substitutions at a freq<4 as 'NoEdit'---rather than discard artifacts that do not contain Indel and are observed at a >4% freq
#
#         #####_____________________________________________________________________________######
# #
# #
#     # print(data[keys])
# #
#     sample_edited = data[keys][data[keys]['sus_artifact'].str.contains('real')].copy()
#
#     artifact_sum = data[keys][data[keys]['sus_artifact'].str.contains('artifact')].copy() #<----------------- Artifact summary
#
#     # data[keys].to_csv(keys+r'_logic.csv', sep=',', header=True,index=False) #<-------- Check suspected artifact logic
#
#    # # __Calculated Variables for New Header____#
#     Proper_reads = int(header_dict[keys][1].split(',')[1])
#     Proper_reads_edit = sum((sample_edited.Hit.astype(int)))
#     Total_Read_Pair = int(header_dict[keys][0].split(',')[1])
#     Total_usable_percentage_edit = (float(Proper_reads_edit)/float(Total_Read_Pair))*100
#     sample_edited.iloc[:,2] = (sample_edited.iloc[:,1].astype(int)/Proper_reads_edit)*100
#     QC_level= sample_edited[~sample_edited['Type'].str.contains('NoEdit')]
#
#     if len(QC_level.index) == 0:
#         Sample_Edit_Level = '0.0'
#     else:
#         Sample_Edit_Level = round(sum(QC_level.Percentage),4)
#
#         # sample_level.to_csv(keys+r'_edit_inq.csv', sep=',', header=True,index=False) #<-------- Check edited %
#
#     hit_dict = dict(zip(sample_edited['CIGAR'], sample_edited['Hit'])) #<-------- Creates dictionary of read number associated with CIGAR string
#
#     # #___Create & Append New Header____#
#
#     new_head = pd.DataFrame({'mask': [keys,Total_Read_Pair, Proper_reads_edit,Total_usable_percentage_edit, str(Sample_Edit_Level)],
#                              'name': ['Sample Name:','Total Read Pairs:', 'Proper Read Pairs:','Total usable percentage:','Sample Edit Level:']})
#
#     cols = new_head.columns.tolist()
#     cols = cols[-1:] + cols[:-1] # <--- Reverses the order of columns in DF
#     rehead = new_head[cols]
#
#     summary = summary.append(rehead)
#     print(rehead)
#
#         # #______CIGAR String Block______________________# #
#     changes = QC_level['CIGAR'].tolist()
#     # changes =(changes[26:28]) #<--------------------------- Subset
#     for cig in changes:
#         CIGAR_edit = Cigar(cig)
#         CIG_list = list(CIGAR_edit.items())
#
#         for index, tup in enumerate(CIG_list):
#             ch = CIG_list[index][1]
#             if ch in String_list:
#
#                 slice_l = CIG_list[:index + 1]
#                 slice_m = CIG_list[:index]
#
#                 pos = str(sum([t[0] for t in slice_l]))
#                 pos2 = str(sum([t[0] for t in slice_m])+1)
#                 indelicous = [t[0] for t in slice_m]
#                 variant = slice_l[-1][1]
#                 var_len = slice_l[-1][0]
#
#                 if variant == 'I' or 'D':
#                     pos_list.append(pos2)
#                     varlen_list.append(var_len)
#                 else:
#                     pos_list.append(pos)
#                     varlen_list.append(var_len)
#
#
#                 var_list.append(variant)
#                 cig_pack.append(cig)
#
#                 # print(cig,slice_m,variant,indelcious)
#             else:
#                 pass
#
#     cig_df['CIGAR'] = cig_pack
#     cig_df['Hit'] = cig_df['CIGAR'].map(hit_dict)
#     cig_df['Percentage'] = ((cig_df['Hit'].astype(int)/Proper_reads_edit)*100)
#     cig_df['Position'] = pos_list
#     cig_df['Variant'] = var_list
#     cig_df['Sample_Name'] = str(keys)
#     cig_df['Var_len'] = varlen_list
#
#     cig_alt = cig_df.loc[(cig_df['Variant'] != 'X') & ~cig_df.Position.isin(artifact_pos)] ##<--------- Add Noise "Yes" parameter?
#     cig_alt = cig_alt.drop_duplicates(keep='first')
#
#     cig_alt = cig_alt.sort_values(['Var_len'], ascending=[True])
#
#
#     # Ptable = cig_df.pivot_table(index=['Sample_Name','CIGAR','Variant'], values= ['Percentage'],aggfunc={'Percentage':'first'})
#
#     if Noise == 'Yes':    ## <------------------------------------------ #Necessary to remove edit genotypes with multiple deletions/insertions for composition fig.
#         Ptable_indel = cig_alt
#         Ptable_indel['Variant'] = Ptable_indel.Variant.replace('DD','D')
#         Ptable_indel['Variant'] = Ptable_indel.Variant.replace('II','I')
#         Ptable = Ptable_indel.pivot_table(index=['Sample_Name','CIGAR','Variant','Var_len'], values= ['Percentage'],aggfunc={'Percentage':'first'})
#     else:
#
#         Ptable = cig_df.pivot_table(index=['Sample_Name','CIGAR','Variant','Var_len'], values= ['Percentage'],aggfunc={'Percentage':'first'})
#
#     # if round(Sample_Edit_Level,0) == 0.0:
#     if Sample_Edit_Level == '0.0':
#
#         d = {'Sample_Name': [keys], 'Variant': ['NoEdit'],'Percentage':[100]}
#         No_Edit_DF = pd.DataFrame(data=d)
#         pivot_dict[keys] = No_Edit_DF
#
#
#     else:
#
#         pivot_summary = pivot_summary.append(Ptable)
#         Ptable2 = pd.DataFrame(pivot_summary.to_records())
#         Ptable3 = Ptable2.groupby(['Percentage','Sample_Name','CIGAR'])['Variant'].apply(lambda x: "{%s}" % ', '.join(x)) #<----- combines variants classes into one string
#         Ptable4 = pd.DataFrame(Ptable3)
#         Ptable4 = pd.DataFrame(Ptable4.to_records())
#
#         if Noise == 'Yes':
#             Ptable4.Variant = Ptable4.Variant.replace('{D, X}','{D}')
#         else:
#             pass
#
#         Ptable5 = Ptable4.pivot_table(index=['Sample_Name','Variant'], values= ['Percentage'],aggfunc={'Percentage':'sum'})
#         Ptable5 = pd.DataFrame(Ptable5.to_records())
#         Ptable5['Variant'] = Ptable5['Variant'].str.replace('\W', '')
#
#         NoEdit_table = Ptable5.pivot_table(index=['Sample_Name'], values= ['Percentage'],aggfunc={'Percentage':'sum'})
#         NoEdit_table = pd.DataFrame(NoEdit_table.to_records())
#         NoEdit_table['Percentage'] = 100 - NoEdit_table['Percentage'].astype(float)
#         NoEdit_table['Variant'] = 'NoEdit'
#         NoEdit_table = NoEdit_table[['Sample_Name','Variant','Percentage']]
#
# # # #   #___________Export Block_______________________________#
# #     artifact_sum.to_csv(keys+r'_artifact_sum.csv', sep=',', header=True,index=False)
# # #
# #     sample_edited.to_csv(keys+r'_QC_Edited_pre.csv', sep=',', header=True,index=False)
# #     rehead.to_csv('rehead.csv', sep=',', header=False,index=False)
# #
# #     with open(keys+r'_QC.csv', 'w') as output:
# #       for f in ['rehead.csv',keys+r'_QC_Edited_pre.csv']:
# #         output.write(''.join([line for line in open(f).readlines() if line.strip()]))
# #
# #
# #     cig_alt.to_csv(keys+r'_CIGAR_summary.csv', sep=',', header=True,index=False)
# # #
# # # # #_____Summary Block_______________________________#
# # summary.to_csv(Locus+r'_summary_knockout.csv', sep=',', header=False,index=False)
#
#
# try: No_Edit_DF
# except NameError: No_Edit_DF = None
#
#
# if No_Edit_DF is None:
#     pivot_dict[keys] = pd.concat([Ptable5,NoEdit_table],axis=0)
#     CompSum = pd.concat(pivot_dict.values(), ignore_index=True)
#     CompSum.to_csv(Locus+r'_pivot_summary.csv', sep=',', header=True,index=False)
#     print("Edit Present")
#
# else:
#     pivot_dict[keys] = pd.concat([No_Edit_DF,Ptable5,NoEdit_table],axis=0)
#     CompSum = pd.concat(pivot_dict.values(), ignore_index=True)
#     CompSum.to_csv(Locus+r'_pivot_summary.csv', sep=',', header=True,index=False)
#     print("No Edit Present")
