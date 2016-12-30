#!usr/bin/env python
# -*- coding: utf-8 -*-

#Load modules
import re
import sys
import glob
from collections import namedtuple
from matplotlib_venn import venn2
import matplotlib.pyplot as plt

directory = '/media/benflies/KING_BEN1/varscan'

vcf_file_list = [f for f in glob.iglob(directory+"/*.vcf")]

set_tst = []
set_hpx = []

for VCF_FILE in vcf_file_list:

    basename = re.sub('.vcf$','',VCF_FILE)
    basename = re.sub(directory,'',basename)
    basename = re.sub('_tst','',basename)
    basename = re.sub('_hpx','',basename)

    IntervalColumns = namedtuple('vcf', ['chr', 'pos', 'ref', 'alt', 'sample'])
    intervals_list = []

    if VCF_FILE:
        with open(VCF_FILE, "r") as fin:
            for line in fin.readlines():
                line = line.rstrip('\n')
                line = re.split(r'\t+', line.rstrip('\t'))

                if len(line) == 10 or len(line)==8 or len(line)==11:
                    if line[0] != '#CHROM':
                        if line[3] == 'A' or line[3] == 'T' or line[3] == 'C' or line[3] == 'G':
                            if line[4] == 'A' or line[4] == 'T' or line[4] == 'C' or line[4] == 'G':

                                line[7] = line[7].split(';')
                                #vcf_line = IntervalColumns(*(line[0], line[1], line[3], line[4], basename))
                                #intervals_list.append(vcf_line)
                                variant_ID = basename+'_'+line[0]+'_'+line[1]+'_'+line[3]+'_'+line[4]
                                #print variant_ID

                                if 'tst' in VCF_FILE:
                                    set_tst.append(variant_ID)
                                if 'hpx' in VCF_FILE:
                                    set_hpx.append(variant_ID)


print set_tst
print('#####')
print set_hpx



set_tst = set(set_tst)
set_hpx = set(set_hpx)

venn2([set_tst, set_hpx], ('TST15', 'Haloplex'))
plt.show()
