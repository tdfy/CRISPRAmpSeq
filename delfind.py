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

bam = sys.argv[1]
path = sys.argv[2]

CIGAR_list = []
String_list = [1,2]
pos_list = []
varlen_list = []
var_list = []
cig_pack = []
ch_list = []


cig_df = pd.DataFrame()


bamFP = pysam.Samfile(path + bam, "rb");

for read in bamFP:
    if( not( read.is_unmapped ) ):
        cig_string = read.cigar
        CIGAR_list.append(cig_string)

        for index, tup in enumerate(cig_string):
            ch = cig_string[index][0]
            if ch in String_list:

                slice_l = cig_string[:index + 1]
                slice_m = cig_string[:index]

                # pos = str(sum([t[1] for t in slice_l]))
                pos2 = str(sum([t[1] for t in slice_m]))
                # indelicous = [t[0] for t in slice_m]
                var_len = slice_l[-1][1]
                variant = slice_l[-1][0]

                # print(cig_string,"slice_m------>",slice_m,pos2)

                if variant == 1 or 2:
                    pos_list.append(pos2)
                    varlen_list.append(var_len)
                    cig_pack.append(cig_string)
                    ch_list.append(ch)
                else:
                    pass


                var_list.append(variant)



    else:
        pass

cig_df['CIGAR'] = cig_pack
cig_df['Position'] = pos_list
cig_df['Variant'] = var_list
cig_df['Var_len'] = varlen_list
cig_df['ch'] = ch_list

cig_df['Variant'] = cig_df['Variant'].map({1:'I', 2:'D'})

# print(cig_df)
insertion = cig_df.loc[(cig_df['Variant'] == 'I')]
Ptable_ins = insertion.pivot_table(index=['Variant','Var_len'], values= ['ch'],aggfunc={'count'})
Ptable2_ins = pd.DataFrame(Ptable_ins)
Ptable3_ins = pd.DataFrame(Ptable2_ins.to_records())


deletion = cig_df.loc[(cig_df['Variant'] == 'D')]
Ptable_del = deletion.pivot_table(index=['Variant','Var_len'], values= ['ch'],aggfunc={'count'})
Ptable2_del = pd.DataFrame(Ptable_del)
Ptable3_del = pd.DataFrame(Ptable2_del.to_records())
Ptable3_del["('ch', 'count')"] = 0 - Ptable3_del["('ch', 'count')"]

df_cat = pd.concat([Ptable3_ins,Ptable3_del],0)


#_____Export Block________________#

df_cat.to_csv(os.path.join(path+bam+r'.csv'), sep=',', header=True,index=True)

# cig_df = cig_df.iloc[0:5,]
cig_df['Hit'] = 1/len(CIGAR_list)*100
cig_df = cig_df[['Position','Variant','Var_len','ch','Hit']]
cig_df.to_csv(os.path.join(path+bam+r'CIGAR_STRING.csv'), sep=',', header=True,index=True)
#
# # pink = cig_df.loc[(cig_df['Var_len'] == 1) & (cig_df['Variant'] == 'D')]
#
# # pink.to_csv(os.path.join(path+bam+r'CIGAR_STRING.csv'), sep=',', header=True,index=True)
#
print("sample", bam)
print("Read Count",len(CIGAR_list))
print("Script is Finished...")
sys.exit()
