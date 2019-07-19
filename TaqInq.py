


import sys
import numpy as np
import pandas as pd
import glob
import os
from cigar import Cigar
import collections
pd.set_option('display.width', 1000)

path = sys.argv[1]
os.chdir(path)

filelist = glob.glob(path+r'/chr14-105062883.Whole_Blood.D21*.ff')

name_list = []
FASTA_list = []

print(filelist)
#
# file = open(filelist[0],"w")
# for x in filelist[0]:
#     if x.startswith('<'):
#         file.write(x)
#
# print(file)

for i in filelist:
    i = pd.read_csv(i, sep=" ", dtype =object)
    #
    # for row in i:
    #     if row.startswith('<'):
    #         print(row)
    #         # row.append(name_list)
    #     else:
    #         # row.append(FASTA_list)
    #         pass
    #

    print(i)
