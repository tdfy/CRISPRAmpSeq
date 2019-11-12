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


# #Modules/Libraries:--------------------------------------------------------------------------------------------------------------

import sys
import numpy as np
import pandas as pd
from cigar import Cigar
import collections
import argparse
from tabulate import tabulate
from collections import defaultdict

# # Global Declarations:--------------------------------------------------------------------------------------------------------------
data = {}
header_dict = {}
pivot_dict = {}
pivot_summary = pd.DataFrame()
String_list = ['I','D','X']
B2M_SNPlist = []
New_File_Name = []
New_File_Dir = []
structure = {}
control_dict = {}
ctrl_df_dict = {}
KO_cat = []
rehead_dict = {}
cig_df_dict = {}
compress_dict = {}
CIGAR_sum_dict = {}

germSNP= defaultdict(lambda: 'NoTarg')
germSNP = {'TRAC':['198=1X36=','1S198=1X36=','197=1X36=','198=1X35=1X1S','188=1X46='], 'CIITA':['52=1X176=2S','1X51=11X176=2S','1X51=1X176=2S'],
'PIK':['39=1D182='],'B2M':['NA']}

# # User Variables:--------------------------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Application to call CRISPR edits from amplicon sequencing libraries",add_help= False)
parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,help=print("""
Configuration Table Example:
""",
tabulate([['/mnt/c/Export/Amp_test/PACE_21','90930_B2M_1.output.csv','B2M','1','T','PACE_20_B','1','1','Y'], ['/mnt/c/Export/Amp_test/PACE_21','90930_B2M_2.output.csv','B2M','2','C','PACE_20_B','2','2','Y']],
headers=['Directory', 'File_Name', 'Target','Sample', 'Design', 'Experiment_Name','Time_Point','Chron_Order','Noise'])))

parser.add_argument('-p','--path', help='path to configuration file...', required = True)
parser.add_argument('-f','--file', help='name of configuration file.[CSV]', required = True)

args = vars(parser.parse_args())

path = args['path']
config = args['file']


configed = pd.read_csv(path +'/'+ config, sep=",", dtype =object)
###-----------------------------------------------------------------------------------------------------------------------------###

for line in configed['File_Name']:
    if line.endswith("v") == True:
        line = line[:-11]
        New_File_Name.append(line)
    else:
        New_File_Name.append(line)

configed['File_Name'] = New_File_Name


####===============================================###

for line in configed['Directory']:
    if line.endswith("/") == False:
        line = str(line) + '/'
        New_File_Dir.append(line)
    else:
        New_File_Dir.append(line)

configed['Directory'] = New_File_Dir

####===============================================###

Targs = configed['Target'].unique()

configed['Full_Name'] = configed['Experiment_Name'] +'_'+ configed['Time_Point']

sample_dict = dict(zip(configed['Sample'], configed['Full_Name']))


for Tag, comp in configed.groupby('Target'):
    structure[Tag] = comp.sort_values(['Design'], ascending=[False])

    control_df = comp.loc[(comp['Design'] == 'C')]
    control_dict[Tag] = control_df['File_Name'].tolist()

    dir = comp['Directory'].iloc[0]

    for samp in structure[Tag]['File_Name']:
        data[samp] = pd.read_csv(dir + samp + '.output.csv', sep=",", dtype =object,skiprows=4)

        ###------ --------------------------------------------------------------------------------------------------------###  ^^^^ Opens output csv's as df's & generates list of controls per target group
        if len(control_dict[Tag]) == 1:
            data[samp] = data[samp][data[samp]['Percentage'].astype(float) >= .01] #<--------- threshold logic (.01)
        else:
            data[samp] = data[samp][data[samp]['Percentage'].astype(float) >= .1] #<--------- threshold logic

        ###------ --------------------------------------------------------------------------------------------------------###  ^^^^ filters reads by proportion criteria based on number of controls

        data[samp].iloc[:,4] = np.where((data[samp].iloc[:,5].isin(germSNP[comp['Target'].iloc[0]])), 'NoEdit', data[samp].iloc[:,4])  ##<------- labels 'Edit' column according to list of germline SNP CIGAR strings; if present in germSNP dict reassigned as 'NoEdit'

         ###---------------------------------------------------------------------------------------------------------------###  ^^^^ Removes known germline variants {convert to dictionary}

for file in control_dict:
    CIGAR_cat = []
    ctrl = control_dict[file]

    if  len(ctrl) > 1:
        for x in ctrl:
            control_samp = data[x].iloc[:,5]
            CIGAR_frame = control_samp.to_frame()     ### converts series to df in order to concat into one DF
            CIGAR_cat.append( CIGAR_frame )     ### appends sub-dfs
            ctrl_df_dict[file] = pd.concat(CIGAR_cat,axis=1,)     ### concats to one DF
    else:
        x = ctrl[0]
        ctrl_df_dict[file] = data[x].iloc[:,5].to_frame()

    ###--------------------------------------------------------------------------------------------------------------------------###


    # ___Call Intersections___#

for samp in structure:
    filelist = structure[samp]['File_Name']

    if len(ctrl_df_dict[samp].columns) == 2:
        u = set.intersection(set(ctrl_df_dict[samp].iloc[:,0]),set(ctrl_df_dict[samp].iloc[:,1]))
        print("Number of Controls:",len(ctrl_df_dict[samp].columns))

    elif len(ctrl_df_dict[samp].columns) == 3:
            u1 = set.intersection(set(ctrl_df_dict[samp].iloc[:,0]),set(ctrl_df_dict[samp].iloc[:,1]))
            u2 = set.intersection(set(ctrl_df_dict[samp].iloc[:,0]),set(ctrl_df_dict[samp].iloc[:,2]))
            u3 = set.intersection(set(ctrl_df_dict[samp].iloc[:,1]),set(ctrl_df_dict[samp].iloc[:,2]))
            print("Number of Controls:",len(ctrl_df_dict[samp].columns))

            u= u1|u2|u3
            print("Length of Artifacts:",len(u1),len(u2),len(u3))

    elif len(ctrl_df_dict[samp].columns) == 1:
            u = []
            cig_df = pd.DataFrame()
            pos_list = []
            var_list = []
            cig_pack = []

            changes = ctrl_df_dict[samp].iloc[:,0].tolist()
            for cig in changes:
                CIGAR_edit = Cigar(cig)
                CIG_list = list(CIGAR_edit.items())

                for index, tup in enumerate(CIG_list):
                    ch = CIG_list[index][1]
                    if ch == 'X':

                        slice_l = CIG_list[:index + 1]
                        pos = str(sum([t[0] for t in slice_l]))
                        variant = slice_l[-1][1]
                        pos_list.append(pos)

                    else:
                        pass

            artifact_pos = [item for item, count in collections.Counter(pos_list).items() if count > 50] #<--------- Threshold for substitution pileup

            print(samp,artifact_pos)

            new_list = [int(x)-1 for x in artifact_pos]

            for file in filelist:
                data[file]['sus_artifact'] = 'real'
                for row in file:
                    for position in new_list:
                        data[file]['Type'] = np.where((data[file]['Type']=='S') & (position in data[file]['CIGAR']),'NoEdit',data[file]['Type'])


    else:
        print("No Control Present")
        exit()

    #### _______Export_list________________________#

    with open(dir +comp['Target'].iloc[0]+r"_artifacts.txt", "w") as output:
        output.write(str(u))


    #####__________________________QC Analysis_____________________________________#

    for file in filelist:

        cig_df = pd.DataFrame()
        pos_list = []
        var_list = []
        cig_pack = []
        varlen_list = []

        with open(dir + file + '.output.csv' ) as myfile:
            header_dict[file] = [next(myfile) for x in range(4)]

        QC = pd.read_csv(dir + file + '.output.csv' , sep=",", dtype =object,skiprows=4)

        if len(u) == 1:
                data[file]['sus_artifact'] = 'real'
        else:

            data[file]['sus_artifact'] = np.where((data[file].iloc[:,5].isin(u)) & (data[file].iloc[:,4] != 'NoEdit')  |
                (data[file].iloc[:,5] == '*') ,'artifact','real')


            data[file]['sus_artifact'] = np.where((data[file].iloc[:,5].str.contains('D|I')==False) & (data[file].iloc[:,2].astype(float) > 4) &
                (data[file].iloc[:,4] != 'NoEdit'),'artifact',data[file]['sus_artifact']) #<----------------- Removes substitution artifacts that do not contain Indel and are observed at a >4% freq


            data[file].iloc[:,4] = np.where((data[file].iloc[:,5].isin(u)) & (data[file].iloc[:,4] != 'NoEdit')  |
                (data[file].iloc[:,5] == '*') ,'NoEdit',data[file].iloc[:,4])

            # data[file].iloc[:,4] = np.where((data[file].iloc[:,5].str.contains('D|I')==False),'NoEdit', data[file].iloc[:,4]) #<----------------- Redefines reads that are substitutions at a freq<4 as 'NoEdit'---rather than discard artifacts that do not contain Indel and are observed at a >4% freq


            #####_____________________________________________________________________________######


        sample_edited = data[file][data[file]['sus_artifact'].str.contains('real')].copy()

        artifact_sum = data[file][data[file]['sus_artifact'].str.contains('artifact')].copy() #<----------------- Artifact summary

       # # __Calculated Variables for New Header____#
        Proper_reads = int(header_dict[file][1].split(',')[1])
        Proper_reads_edit = sum((sample_edited.Hit.astype(int)))
        Total_Read_Pair = int(header_dict[file][0].split(',')[1])
        Total_usable_percentage_edit = (float(Proper_reads_edit)/float(Total_Read_Pair))*100
        sample_edited.iloc[:,2] = (sample_edited.iloc[:,1].astype(int)/Proper_reads_edit)*100
        QC_level= sample_edited[~sample_edited['Type'].str.contains('NoEdit')]


        if len(QC_level.index) == 0:
            Sample_Edit_Level = '0.0'
        else:
            Sample_Edit_Level = round(sum(QC_level.Percentage),4)

        hit_dict = dict(zip(sample_edited['CIGAR'], sample_edited['Hit']))

        # #___Create & Append New Header____#

        new_head = pd.DataFrame({'mask': [file,Total_Read_Pair, Proper_reads_edit,Total_usable_percentage_edit, str(Sample_Edit_Level)],
                                 'name': ['Sample Name:','Total Read Pairs:', 'Proper Read Pairs:','Total usable percentage:','Sample Edit Level:']})

        cols = new_head.columns.tolist()
        cols = cols[-1:] + cols[:-1] # <--- Reverses the order of columns in DF
        rehead = new_head[cols]


        rehead_dict[file] = rehead

        # print(rehead)

            # #______CIGAR String Block______________________# #
        changes = QC_level['CIGAR'].tolist()
        for cig in changes:
            CIGAR_edit = Cigar(cig)
            CIG_list = list(CIGAR_edit.items())

            for index, tup in enumerate(CIG_list):
                ch = CIG_list[index][1]
                if ch in String_list:

                    slice_l = CIG_list[:index + 1]
                    slice_m = CIG_list[:index]

                    pos = str(sum([t[0] for t in slice_l]))
                    pos2 = str(sum([t[0] for t in slice_m])+1)
                    variant = slice_l[-1][1]
                    var_len = slice_l[-1][0]

                    if variant == 'I' or 'D':
                        pos_list.append(pos2)
                        varlen_list.append(var_len)
                    else:
                        pos_list.append(pos)
                        varlen_list.append(var_len)


                    var_list.append(variant)
                    cig_pack.append(cig)

                else:
                    pass

        cig_df['CIGAR'] = cig_pack
        cig_df['Hit'] = cig_df['CIGAR'].map(hit_dict)
        cig_df['Percentage'] = ((cig_df['Hit'].astype(int)/Proper_reads_edit)*100)
        cig_df['Position'] = pos_list
        cig_df['Variant'] = var_list
        cig_df['Sample_Name'] = str(file)
        cig_df['Var_len'] = varlen_list


        if (len(u) == 1) & (comp['Noise'].iloc[0] == 'Yes'):
            cig_alt = cig_df.loc[(cig_df['Variant'] != 'X') & ~cig_df.Position.isin(artifact_pos)]
            cig_alt = cig_alt.drop_duplicates(keep='first')

            cig_alt = cig_alt.sort_values(['Var_len'], ascending=[True])

        else:
            cig_alt = cig_df

        if comp['Noise'].iloc[0] == 'Yes':    ## <------------------------------------------ #Necessary to remove edit genotypes with multiple deletions/insertions for composition fig.
            Ptable_indel = cig_alt
            Ptable_indel['Variant'] = Ptable_indel.Variant.replace('DD','D')
            Ptable_indel['Variant'] = Ptable_indel.Variant.replace('II','I')
            prePtable = Ptable_indel.pivot_table(index=['Sample_Name','CIGAR','Variant','Var_len'], values= ['Percentage'],aggfunc={'Percentage':'first'})

            if prePtable.empty:
                pass
            else:
                Ptable = pd.DataFrame(prePtable.to_records())
        else:

            prePtable = cig_df.pivot_table(index=['Sample_Name','CIGAR','Variant','Var_len'], values= ['Percentage'],aggfunc={'Percentage':'first'})

            if prePtable.empty:
                pass
            else:
                Ptable = pd.DataFrame(prePtable.to_records())


        if prePtable.empty:
                d = {'Sample_Name': [file], 'CIGAR':'NA', 'Variant': ['NoEdit'],'Var_len':0,'Percentage':[100]}
                Ptable = pd.DataFrame(data=d)
        else:
            pass

        cig_df_dict[file] = Ptable

        ##------------- New Export Block --------------------##

        cig_alt.to_csv(dir+file+r'_CIGAR_summary.csv', sep=',', header=True,index=False)

        sample_edited.to_csv(dir+file+r'_QC_Edited_pre.csv', sep=',', header=True,index=False)
        rehead.to_csv(dir+'rehead.csv', sep=',', header=False,index=False)

        with open(dir+file+r'_QC.csv', 'w') as output:
          for f in [dir + 'rehead.csv',dir + file+r'_QC_Edited_pre.csv']:
            output.write(''.join([line for line in open(f).readlines() if line.strip()]))

            ##-----------------------------------------------##

for loci in Targs:
    key_list =  list(cig_df_dict.keys())
    sorted_keys = [k for k in key_list if loci in k]    ####<<<<------------ creates sub set list by target loci
    cig_select = {akey:cig_df_dict[akey] for akey in sorted_keys if akey in cig_df_dict} ####<<<<------------ creates sub set dictionary from cig_df_dict based on sorted_keys list
    Ptable2 = pd.concat(cig_select.values(),ignore_index=True)
    Ptable3 = Ptable2.groupby(['Percentage','Sample_Name','CIGAR'])['Variant'].apply(lambda x: "{%s}" % ', '.join(x)) #<----- combines variants classes into one string
    Ptable4 = pd.DataFrame(Ptable3)
    Ptable4 = pd.DataFrame(Ptable4.to_records())
    Ptable4['Variant'] = Ptable4['Variant'].str.replace('\W', '') #<----- removes comma sep and {}


    if comp['Noise'].iloc[0] == 'Yes':
        Ptable4.Variant = np.where((Ptable4.Variant.str.startswith('D')),'D',Ptable4.Variant)
        Ptable4.Variant = np.where((Ptable4.Variant.str.startswith('I')),'I',Ptable4.Variant)

    else:
        pass

    CIGAR_sum_dict[loci] = Ptable4

    Ptable5 = Ptable4.pivot_table(index=['Sample_Name','Variant'], values= ['Percentage'],aggfunc={'Percentage':'sum'})
    Ptable5 = pd.DataFrame(Ptable5.to_records())

    NoEdit_table = Ptable5.pivot_table(index=['Sample_Name'], values= ['Percentage'],aggfunc={'Percentage':'sum'})
    NoEdit_table = pd.DataFrame(NoEdit_table.to_records())
    NoEdit_table['Percentage'] = 100 - NoEdit_table['Percentage'].astype(float)
    NoEdit_table['Variant'] = 'NoEdit'
    NoEdit_table = NoEdit_table[['Sample_Name','Variant','Percentage']]

    combo = pd.concat([Ptable5, NoEdit_table],axis=0)

    compress_dict[loci]= combo


pivot_sum = pd.concat(compress_dict.values(),ignore_index=True)
pivot_sum = pivot_sum[pivot_sum['Percentage'] != 0 ]

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

    # #_____Summary Block_______________________________#

summary = pd.concat(rehead_dict.values(),ignore_index=True)
summary.to_csv(dir + comp['Experiment_Name'].iloc[0]+r'_summary_knockout.csv', sep=',', header=False,index=False)

# #__________________KO Summary Formating____________________________#

hope = np.array_split(summary, (len(sample_dict)*len(Targs)))

for i in hope:
    i_tran= i.transpose()
    i_tran.columns = ['Sample Name:','Total Read Pairs:','Proper Read Pairs:','Total usable percentage:','Sample Edit Level:']
    KO_cat.append(i_tran) ### appends sub-dfs

KO_df = pd.concat(KO_cat,axis=0)     ### concats to one DF
KO_df= KO_df[KO_df['Sample Name:'].str.contains('_')]#<-------- Assumes all Specimen Names contain '_'

KO_df['Samp_No'] = KO_df['Sample Name:'].str.split('_').str.get(2)
KO_df['Locus'] = KO_df['Sample Name:'].str.split('_').str.get(1)

KO_df = KO_df.sort_values(['Samp_No','Locus'], ascending=[True,True])

KO_df["Insertion/Deletion_Level_(%)"] = KO_df['Sample Name:'].map(indel_dict)

KO_df["Specimen_Name"] = KO_df['Samp_No'].map(sample_dict)

KO_df = KO_df[['Specimen_Name','Locus','Total Read Pairs:','Proper Read Pairs:','Total usable percentage:','Sample Edit Level:',"Insertion/Deletion_Level_(%)"]]

KO_df['Total usable percentage:'] = round(KO_df['Total usable percentage:'].astype(float),2)
KO_df['Sample Edit Level:'] = round(KO_df['Sample Edit Level:'].astype(float),2)

# #_______________Output Block_____________________________________________# #
KO_df.to_csv(dir + comp['Experiment_Name'].iloc[0]+'.csv', sep=',', header=True,index=False)

print(KO_df)
