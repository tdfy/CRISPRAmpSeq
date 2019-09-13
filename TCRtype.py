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
import pybedtools

bam = sys.argv[1]
path = sys.argv[2]

CIGAR_list = []
String_list = [1,2]
pos_list = []
varlen_list = []
var_list = []
cig_pack = []
# ch_list = []
CIG_dict = {}
ID = []
start_coo = []
seq_list = []
read_list = []

cig_df = pd.DataFrame()

trans_dict = {}

genecode = 'trans_annotate.csv'
annotate = pd.read_csv(path +'/'+ genecode, sep=",", dtype =object)
anno_dict = dict(zip(annotate['name'], annotate['Gene']))

annotate['chr'] = "7"
annotate_bed = pybedtools.BedTool.from_dataframe(annotate[['chr','txStart','txEnd','name']])

# ENST00000436911 = pd.read_csv(path + '/' +'ENST00000436911.6.bed',sep=',',dtype=object)

bamFP = pysam.Samfile(path + bam, "rb");

for read in bamFP:
    if( not( read.is_unmapped ) ):
        cig_string = read.cigar
        CIGAR_list.append(cig_string)
        cigarette = read.cigarstring
        pos_uni = read.get_reference_positions() # full_length=True
        cig_pos = pos_uni[0]
        seq = read.query_alignment_sequence
        name = read.query_name


        for index, tup in enumerate(cig_string):
            if cig_string[0][0] == 4 or cig_string[0][0] == 5 :
                del cig_string[0]
            else:
                pass

            ch = cig_string[index][0]

            if ch in String_list:

                slice_l = cig_string[:index + 1]
                slice_m = cig_string[:index]

                pos2 = str(sum([t[1] for t in slice_m]))

                var_len = slice_l[-1][1]
                variant = slice_l[-1][0]

                if variant == 1 or 2:
                    pos_list.append(pos2)
                    varlen_list.append(var_len)
                    cig_pack.append(cig_string)
                    start_coo.append(cig_pos)
                    seq_list.append(seq)
                    read_list.append(name)
                    var_list.append(variant)
                    ID.append(cigarette)

                    # print(cigarette,cig_pos,len(pos_uni),pos2,seq)

                    # ch_list.append(ch)
                else:
                    pass


                # var_list.append(variant)
                # ID.append(cigarette)



    else:
        pass


index_list = list(range(1, len(CIGAR_list)))
cig_df['CIGAR'] = cig_pack

cig_df['Var_len'] = varlen_list

cig_df['Position'] = pos_list
cig_df['Read_Start'] = start_coo
cig_df['Read_Start'] = cig_df['Read_Start'].astype(int) + 38240024
cig_df['Position_Start'] = cig_df['Position'].astype(int) + cig_df['Read_Start'].astype(int) #+ 38240024 #'164419787'  #<------------------- make starting reference coordinate a variable "mmPIGA"
cig_df['Position_Stop'] = cig_df['Position_Start'].astype(int) + cig_df['Var_len'].astype(int)

cig_df['Variant'] = var_list
# cig_df['ch'] = ch_list

cig_df['Variant'] = cig_df['Variant'].map({1:'I', 2:'D'})
cig_df['ID'] = ID
cig_df['FASTA'] = seq_list
cig_df['Read_ID'] = read_list

# cig_df['Clone_ID'] = cig_df.Read_Start.astype(str) + cig_df.Var_len.astype(str)+ '_' + cig_df.Position.astype(str) + cig_df.Variant
cig_df['Clone_ID'] = cig_df.Position_Start.astype(str)+ '_' + cig_df.Var_len.astype(str)+ '_' + cig_df.Variant #cig_df.Position.astype(str)



cig_df['chr'] = '7'

samp = cig_df[['chr','Position_Start','Position_Stop','Clone_ID']]


####____________BedTools Intersection______________________________________________####
samp_bed = pybedtools.BedTool.from_dataframe(samp)
N = pybedtools.BedTool.intersect(samp_bed,annotate_bed,wb=True)

N = N.to_dataframe()

trans_dict = dict(zip(N['name'], N['thickEnd']))

cig_df['Transcript'] = cig_df['Clone_ID'].map(trans_dict)
cig_df['Transcript'] = cig_df['Transcript'].fillna('intronic')
cig_df['Annotation'] = cig_df['Transcript'].map(anno_dict)
cig_df['Annotation'] = cig_df['Annotation'].fillna('intronic')

cig_df['Clone_Num'] = cig_df.groupby(['Var_len','Position_Start','Variant','Transcript'])['Clone_ID'].transform('count')
cig_df['Clonotype'] = cig_df.groupby(['Var_len','Position_Start','Variant','Transcript'])['Clone_ID'].transform('first')
#

cig_df.to_csv(path+'/'+ bam +'_ANALYSIS.csv', sep=',', header=True,index=False)

select = cig_df[['Transcript','Annotation', 'Clonotype','Clone_Num']]

# summary = select.groupby(['Transcript','Annotation', 'Clonotype'])['Clone_Num'].transform('first')
summary = select.pivot_table(index=['Transcript','Annotation', 'Clonotype'], values= ['Clone_Num'],aggfunc={'first'})
summary = pd.DataFrame(summary.to_records())
summary.to_csv(path+'/'+ bam +'_SUM.csv', sep=',', header=True,index=False)
print(summary)

print("Total Number of Clones",len(cig_df['CIGAR']))

print("Total Number Reads",len(CIGAR_list))

print("Total Number Clonotypes",len(summary['Transcript']))

###==================================================================================####


# for i in annotate['name']:
#     cig_df['Gene'] = np.where((cig_df['Position_Start'] >= annotate['txStart']) & (cig_df['Position_Start'] <= annotate['txStop']), annotate.name,'NaN')
#
#
# print(cig_df)
# #
# print(len(cig_df['Clonotype'].unique()))


# insertion = cig_df.loc[(cig_df['Variant'] == 'I')]
# Ptable_ins = insertion.pivot_table(index=['Variant','Var_len','ID','Position','Position_Start','Position_Stop'], values= ['ch'],aggfunc={'count'})
# Ptable2_ins = pd.DataFrame(Ptable_ins)
# Ptable3_ins = pd.DataFrame(Ptable2_ins.to_records())
#
#
# deletion = cig_df.loc[(cig_df['Variant'] == 'D')]
# Ptable_del = deletion.pivot_table(index=['Variant','Var_len','ID','Position','Position_Start','Position_Stop'], values= ['ch'],aggfunc={'count'})
# Ptable2_del = pd.DataFrame(Ptable_del)
# Ptable3_del = pd.DataFrame(Ptable2_del.to_records())
# Ptable3_del["('ch', 'count')"] = 0 - Ptable3_del["('ch', 'count')"]
#
# df_cat = pd.concat([Ptable3_ins,Ptable3_del],0)
#
#
# #_____Export Block________________#
#
# df_cat.to_csv(os.path.join(path+bam+r'.csv'), sep=',', header=True,index=True)
#
# # cig_df = cig_df.iloc[0:5,]
# cig_df['Hit'] = 1/len(CIGAR_list)*100
# cig_df = cig_df[['Position','Position_Start','Position_Stop','Variant','Var_len','ch','Hit','ID']]
# cig_df.to_csv(os.path.join(path+bam+r'CIGAR_STRING.csv'), sep=',', header=True,index=True)
# #
# print("sample", bam)
# print("Read Count",len(CIGAR_list))
# print("Total Edits", len(np.unique(cig_df['ID'])))
# print("Script is Finished...")
# sys.exit()
