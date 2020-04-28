#!/usr/bin/python

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
# Author:	Todd D. Yoder
#
# Date:	April 2019
#


# #Modules:--------------------------------------------------------------------------------------------------------------

import sys
import numpy as np
import pandas as pd
import glob
import os
from cigar import Cigar
import collections
import pysam

# bam = sys.argv[1]
path = sys.argv[1]
Target = int(sys.argv[2])
# ref_start = int(sys.argv[3])

bam_list = []
final_sum_df = pd.DataFrame()
sum_dict = {}
String_list = ['I','D']


with open(path+'/'+'bam_list', 'r') as f:
    x = f.readlines()
    for j in x:
        j = j.strip('\n')
        bam_list.append(j)

# print(bam_list)
for bam in bam_list:

    CIGAR_list = []
    # String_list = [1,2]
    pos_list = []
    varlen_list = []
    var_list = []
    cig_pack = []
    ch_list = []
    CIG_dict = {}
    COO_dict = {}
    read_len_dict = {}
    ID = []
    cig_df = pd.DataFrame()
    len_list = []
    query_list = []
    qscore_list = []
    filter_df = pd.DataFrame()
    C_list = []
    indel_list = []
    ref_st_list = []
    coo_list = []
    cig_len = []

    bamFP = pysam.Samfile(path + bam, "rb");

    for read in bamFP:

        if( not( read.is_unmapped ) ):
            cig_string = read.cigar
            CIGAR_list.append(cig_string)
            cigarette = read.cigarstring
            C_list.append(cigarette)
            # print(cigarette)

            # read_len = int(read.infer_read_length())
            read_len = read.reference_length

            len_list.append(read_len)


            query = read.qname
            query_list.append(query)

            qs = read.mapping_quality
            qscore_list.append(qs)

            start = read.reference_start
            ref_st_list.append(start)

    filter_df['Read_Name'] = query_list
    filter_df['Score'] = qscore_list
    filter_df['Score2'] = qscore_list
    filter_df['CIGAR_str'] = CIGAR_list
    filter_df['CIGAR_list'] = C_list
    filter_df['read_len'] = len_list
    filter_df['start_coo'] = ref_st_list
    filter_df['loc'] = filter_df.index

    # CIG_dict = dict(zip(filter_df['Score'], filter_df['Read_Name']))

    # hobby = filter_df.groupby(['Read_Name', 'Score'], as_index=False, sort=False)['Score2'].max()
    # indices = hobby.index.tolist()
    # filter_df = filter_df[filter_df.index.isin(indices)]
    # print(filter_df[filter_df['Read_Name']== 'm54257_200224_201147/32834040/ccs'])

    # hobby = filter_df.groupby(['Read_Name', 'Score'], as_index=False, sort=False)['Score2'].max()

    # c_maxes = filter_df.groupby(['Read_Name']).Score.transform(max)
    filter_df = filter_df.sort_values(['Score'], ascending=[False])

    filter_df = filter_df.drop_duplicates(subset=['Read_Name'],keep='first')
    # filter_df2 = filter_df.loc[filter_df.Score2 == c_maxes]


    CIG_dict = dict(zip(filter_df['CIGAR_list'], filter_df['Read_Name']))
    COO_dict = dict(zip(filter_df['CIGAR_list'], filter_df['start_coo']))
    read_len_dict = dict(zip(filter_df['CIGAR_list'], filter_df['read_len']))



    # print(filter_df)
 ###____Old CIGAR string variant calling___###
    # for c_str in filter_df['CIGAR_str']:
    #     for index, tup in enumerate(c_str):
    #             ch = c_str[index][0]
    #             if ch in String_list:
    #
    #                 slice_l = c_str[:index + 1]
    #                 slice_m = c_str[:index]
    #
    #                 pos2 = str(sum([t[1] for t in slice_m]))
    #
    #                 var_len = slice_l[-1][1]
    #                 variant = slice_l[-1][0]
    #
    #                 if variant == 1 or 2:
    #                     pos_list.append(pos2)
    #                     varlen_list.append(var_len)
    #                     cig_pack.append(cig_string)
    #                     ch_list.append(ch)
    #                 else:
    #                     pass
    #
    #
    #                 var_list.append(variant)
    #                 ID.append(readnamey)
    #
    #     else:
    #         pass
    #

###________________________________________________________##

    for c_str in filter_df['CIGAR_list']:
        CIGAR_edit = Cigar(c_str)
        l = len(CIGAR_edit)
        CIG_list = list(CIGAR_edit.items())

        for index, tup in enumerate(CIG_list):
            if CIG_list[0][1] is 'S':
                CIG_list.remove(CIG_list[0])
            else:
                pass
            ch = CIG_list[index][1]

            if ch in String_list:

                slice_l = CIG_list[:index + 1]
                slice_m = CIG_list[:index]

                pos2 = str(sum([t[0] for t in slice_m])+1)
                variant = slice_l[-1][1]
                var_len = slice_l[-1][0]

                if variant == 'I' or 'D':
                    pos_list.append(pos2)
                    varlen_list.append(var_len)
                    cig_pack.append(CIG_list)
                    ch_list.append(ch)
                    var_list.append(variant)
                    indel_list.append(c_str)
                    cig_len.append(l)
                else:
                    pass


                # var_list.append(variant)
                # ID.append(readnamey)

            else:
                pass
    # index_list = list(range(1, len(CIGAR_list)))

    cig_df['CIGAR'] = cig_pack
    cig_df['Var_len'] = varlen_list

    cig_df['Position'] = pos_list
    # cig_df['Position_Start'] = cig_df['Position'].astype(int) + ref_start  #<------------------- make starting reference coordinate a variable
    # cig_df['Position_Stop'] = cig_df['Position_Start'].astype(int) + cig_df['Var_len'].astype(int)

    cig_df['Variant'] = var_list
    cig_df['ch'] = ch_list

    # cig_df['Variant'] = cig_df['Variant'].map({1:'I', 2:'D'})
    cig_df['ID'] = indel_list
    cig_df['COO_start'] = cig_df['ID'].map(COO_dict)
    cig_df['Read_len'] = cig_df['ID'].map(read_len_dict)
    # cig_df['Read_len'] = cig_len


    cig_df['ID']=cig_df['ID'].map(CIG_dict)   ### <<<<------------------------------------------------------------------------------------------------------------------------------------------------ Added quotes to variable, must replace!!!

    # cig_df['Target'] = Target + ref_start
    cig_df['Target'] = Target -cig_df['COO_start']

    cig_df['Position_Start'] = cig_df['Position'].astype(int)  #<------------------- make starting reference coordinate a variable
    cig_df['Position_Stop'] = cig_df['Position_Start'] + cig_df['Var_len'].astype(int)

    cig_df['Edit_Dis_Start'] = abs(cig_df['Target'] - cig_df['Position_Start'])
    cig_df['Edit_Dis_Stop'] = abs(cig_df['Target'] - cig_df['Position_Stop'])
    cig_df['On_Target'] = np.where((cig_df.Edit_Dis_Start <=10) | (cig_df.Edit_Dis_Stop <=10) | ((cig_df['Target'].astype(int)).between(cig_df['Position_Start'].astype(int), cig_df['Position_Stop'].astype(int))),True,False)
    cig_df['On_Target'] = np.where((cig_df['COO_start'] <=1000) & (cig_df['Read_len'] < Target),False,cig_df['On_Target'])
    cig_df['On_Target'] = np.where(cig_df['COO_start'] > Target, False, cig_df['On_Target'])
    cig_df['Cov_Off'] = np.where((cig_df['COO_start'] <=1000) & (cig_df['Read_len'] < Target),1,0)
    cig_df['Cov_Off'] = np.where(cig_df['COO_start'] > Target, 1, cig_df['Cov_Off'])



    # print(len(np.where((cig_df['COO_start'] <=1000) & (cig_df['Read_len'] < Target),False,cig_df['On_Target'])))


    cig_df.to_csv(os.path.join(path+bam+r'_EXPL.csv'), sep=',', header=True,index=True)
    # print(cig_df)
    insertion = cig_df.loc[(cig_df['Variant'] == 'I')]
    Ptable_ins = insertion.pivot_table(index=['Variant','Var_len','ID','Position','Position_Start','Position_Stop','Target','On_Target'], values= ['ch'],aggfunc={'count'})
    Ptable2_ins = pd.DataFrame(Ptable_ins)
    Ptable3_ins = pd.DataFrame(Ptable2_ins.to_records())

    deletion = cig_df.loc[(cig_df['Variant'] == 'D')]
    Ptable_del = deletion.pivot_table(index=['Variant','Var_len','ID','Position','Position_Start','Position_Stop','Target','On_Target'], values= ['ch'],aggfunc={'count'})
    Ptable2_del = pd.DataFrame(Ptable_del)
    Ptable3_del = pd.DataFrame(Ptable2_del.to_records())
    Ptable3_del["('ch', 'count')"] = 0 - Ptable3_del["('ch', 'count')"]

    df_cat = pd.concat([Ptable3_ins,Ptable3_del],0)

    print(cig_df)

    #_____Export Block________________#

    df_cat.to_csv(os.path.join(path+bam+r'.csv'), sep=',', header=True,index=True)

    # cig_df = cig_df.iloc[0:5,]
    # cig_df['Hit'] = 1/sum(i > Target + 50 for i in len_list)*100
    # cig_df = cig_df[['Position','Position_Start','Position_Stop','Variant','Var_len','ch','ID']]
    cig_df.to_csv(os.path.join(path+bam+r'CIGAR_STRING.csv'), sep=',', header=True,index=True)


    # on_targ = sum(i > Target + 100 for i in len_list)
    off_df = cig_df.loc[cig_df['Cov_Off'] == 1] #number of unique reads off target subtracted from the total reads
    on_targ = len(filter_df )- off_df.ID.nunique()
    print(len(filter_df),on_targ)

    print(bamFP)
    print("sample", bam)
    print("Read Count",len(filter_df))
    # print("Total Edits", len(np.unique(cig_df['ID'])))
    print("Number of reads capturing Target",on_targ)

    df_cat['Count_log'] = df_cat.groupby(['ID'])['On_Target'].transform('unique')
    df_count= df_cat[df_cat['Count_log'].astype(str) != '[False]']

    # tot_edit = df_count.ID.nunique()
    tot_edit =  df_count.ID.nunique()


    max_indel = max(cig_df['Var_len'])
    mean_indel = np.mean(cig_df['Var_len'])
    std_indel = np.std(cig_df['Var_len'])
    mean_read_len = np.mean(len_list)
    first_quant = np.percentile(len_list, 25)
    third_quant = np.percentile(len_list, 75)

    summary = [{"Sample":bam,"Read_Count": len(filter_df),"On_Targ_Cov":on_targ,"On_Targ_Edit":tot_edit,
    "Max_indel":max_indel,"Mean_indel":mean_indel,"Indel_std":std_indel,"Mean_Read_len":mean_read_len,"Len_25":first_quant,"Len_75":third_quant}]
    summary_df = pd.DataFrame(summary)

    sum_dict[bam] = summary_df

    # final_sum_df.append(summary_df)
final =  pd.concat(sum_dict.values(),ignore_index=True)

print(final)

# pd.Series(len_list).plot.hist(grid=True, bins=20, rwidth=0.9,
#                    color='#607c8e')
# plt.title('Commute Times for 1,000 Commuters')
# plt.xlabel('Counts')
# plt.ylabel('Commute Time')
# plt.grid(axis='y', alpha=0.75)

sys.exit()
