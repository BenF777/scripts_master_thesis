#! usr/bin/env python
# -*- coding: utf-8 -*-
"""
This script calculates the percentage of reads above a certain
threshold for each position and plots the result
"""
#import modules
from Bio import SeqIO
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
import pylab as p
import glob
import gzip
import seaborn as sns

### Input

directory_hpx_ffpe = '/home/benflies/NGS_Data/FASTQ/fastq_hpx_ffpe'
directory_hpx_cll = '/home/benflies/NGS_Data/FASTQ/fastq_hpx_cll'
directory_hpx_horizon = '/home/benflies/NGS_Data/FASTQ/fastq_hpx_horizon'

fastq_hpx_ffpe = [f for f in glob.iglob(directory_hpx_ffpe+"/*R1_001.fastq.gz")]
fastq_hpx_cll = [f for f in glob.iglob(directory_hpx_cll+"/*R1_001.fastq.gz")]
fastq_hpx_horizon = [f for f in glob.iglob(directory_hpx_horizon+"/*R1_001.fastq.gz")]

sns.set_context("paper")
plt.figure(figsize=(8, 6))

### Open file & parse over it
for file in fastq_hpx_ffpe:
    print file
    handle = gzip.open(file)
    recs = SeqIO.parse(handle,"fastq")
    d1 = defaultdict(list)
    for rec in recs:
        pos = 0
        for i,qual in enumerate(rec.letter_annotations['phred_quality']):
            d1[pos].append(qual)
            pos = pos + 1
    means = []
    mean = 0
    for key, value in d1.items():
        mean = np.mean(value)
        means.append(mean)
    x = np.arange(len(means))
    plt.plot(x, means, color='red')

for file in fastq_hpx_cll:
    print file
    handle = gzip.open(file)
    recs = SeqIO.parse(handle,"fastq")
    d1 = defaultdict(list)
    for rec in recs:
        pos = 0
        for i,qual in enumerate(rec.letter_annotations['phred_quality']):
            d1[pos].append(qual)
            pos = pos + 1
    means = []
    mean = 0
    for key, value in d1.items():
        mean = np.mean(value)
        means.append(mean)
    x = np.arange(len(means))
    plt.plot(x, means, color='blue')

for file in fastq_hpx_horizon:
    print file
    handle = gzip.open(file)
    recs = SeqIO.parse(handle,"fastq")
    d1 = defaultdict(list)
    for rec in recs:
        pos = 0
        for i,qual in enumerate(rec.letter_annotations['phred_quality']):
            d1[pos].append(qual)
            pos = pos + 1
    means = []
    mean = 0
    for key, value in d1.items():
        mean = np.mean(value)
        means.append(mean)
    x = np.arange(len(means))
    plt.plot(x, means, color='blue')
    
plt.ylim(0, 40)
plt.ylabel('Quality (Phred Score)')
plt.xlabel('Position (bp)')
plt.show()
