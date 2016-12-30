#!usr/bin/env python
# -*- coding: utf-8 -*-

'''
This script plots the normalized coverage distribution per amplicon in a
given directory. The average coverage per target region is normalized by
dividing this value by the overall average coverage of each target sequences.
'''

#Load modules
import pybedtools
import re
import sys
import glob
from collections import namedtuple, defaultdict, OrderedDict
import numpy as np
import matplotlib.pyplot as plt
from operator import itemgetter
import seaborn as sns
import pandas as pd

INTERVALS_BED_A = '/home/benflies/NGS_Data/00100-1407755742_Regions.bed'

directory_A = '/home/benflies/NGS_Data/test'

sns.set_context("paper")
plt.figure(figsize=(8, 6))

fig,ax1 = plt.subplots()
plt.hold = True
boxes = []

bam_file_list_A = [f for f in glob.iglob(directory_A+"/*.bam")]

collected = defaultdict(list)

for file in bam_file_list_A:

    print file

    almnt = pybedtools.BedTool(file)

    IntervalColumns = namedtuple('bed', ['chr', 'start', 'end', 'gene'])
    intervals_list = []

    if INTERVALS_BED_A:
        with open(INTERVALS_BED_A, "r") as fin:
            for line in fin.readlines():
                line = line.rstrip('\n')
                line = re.split(r'\t+', line.rstrip('\t'))

                # For Agilent Haloplex BED files
                if len(line) == 4:
                    line[1] = int(line[1])
                    line[2] = int(line[2])
                    bed_line = IntervalColumns(*(line[0], int(line[1]), int(line[2]), line[3]))
                    intervals_list.append(bed_line)

                # For Illumina Trusight Tumor 15 BED files
                if len(line) == 12:
                    line[1] = int(line[1])
                    line[2] = int(line[2])
                    bed_line = IntervalColumns(*(line[0], int(line[1]), int(line[2]), line[3]))
                    intervals_list.append(bed_line)
    else:
        print("ERROR: Provide an interval list (bed format)")

    coverage_result = almnt.coverage(intervals_list)

    '''
    Print amplicons that have not been amplified at all
    and amplicons that have coverage < 1000x
    '''

    cov_counter = []

    for interval in coverage_result:
        cov_counter.append(float(interval[4]))

    cov_mean = np.mean(cov_counter)

    for interval in coverage_result:
        if interval[3] not in collected.keys():
            collected[interval[3].encode('ascii','ignore')] = {
            'Coverages' : [float(interval[4]) / cov_mean],
            'Median' : 0
            }
        else:
            collected[interval[3].encode('ascii','ignore')]['Coverages'].append(float(interval[4]) / cov_mean)

collected_sorted={}
for key, value in collected.items():
    median = np.median(collected[key]['Coverages'])
    collected[key]['Median'] += median
    collected_sorted = OrderedDict(sorted(collected.items(), key=lambda t: t[1]['Median']))

for key, value in collected_sorted.items():
    boxes.append(collected_sorted[key]['Coverages'])

print collected_sorted

bp = plt.boxplot(boxes,sym='')

xtickNames = plt.setp(ax1, xticklabels = [])

#plt.setp(xtickNames, rotation=90, fontsize=7)
plt.ylabel('Coverage (x)')
plt.xlabel('Target ID')
#plt.title('Comparison of Amplicon Depths Across Samples')
plt.show()
