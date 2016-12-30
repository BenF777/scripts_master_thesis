#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
This script computes run statistics
"""

from __future__ import absolute_import, division, print_function

import sys
import os
import gzip
import csv
import re
import operator
import numpy as n
import pylab as p
import matplotlib.pyplot as plt

from Bio import SeqIO

def get_basename(name, extensions, cut_paired=False):
     """
     Extensions: ['.fastq', '.fastq.gz']
     Cuts of extensions provided and full path,
     If cut_paired, it looks to detect paired signs and cuts those off
     """

     if name is None:
          return None
     else:
          name = os.path.basename(name)
          for ext in extensions:
                name = re.sub("%s$" % ext, '', name)
          if cut_paired is True:
                name = re.sub("_L001_R._001$", '', name)
                name = re.sub("_1$|_2$", '', name)
          return name

directories = sys.argv[1:]

plt.figure()
plt.hold = True
boxes = []

for arg in directories:

    overall_av_counter = 0
    total_reads_counter = 0
    overall_size_counter = 0

    av_list = []
    read_list = []
    size_list = []

    os.chdir("/home")
    path_fastq = arg
    os.chdir(os.path.dirname(os.path.abspath(arg)))
    files = os.listdir(path_fastq)

    #search fastq files:
    f = [os.path.join(root,name)
                 for root, dirs, files in os.walk(path_fastq)
                 for name in files
                 if name.endswith((".fastq", ".fastq.gz"))]

    #collected = {}
    for file in f:

        handle = gzip.open(file)
        reads = 0
        average = []
        av = 0
        sizes = 0

        #count reads & calculate average
        for rec in SeqIO.parse(handle, "fastq"):
            reads += 1
            size = len(rec)
            sizes += len(rec)
            #size_list.append(size)

            if reads % 5000 == 0:
                name = get_basename(file, ['.fastq', '.fastq.gz'])
                print("currently in file %s at read %s" %(name,reads))
            count_qual = 0
            for qual in rec.letter_annotations["phred_quality"]:
                count_qual += qual
            if len(rec) > 0:
                average = count_qual / len(rec)
            else:
                average = 0
            if average > 0:
                av += average
            else:
                av = 0

        if reads > 0:
            if average > 0:
               average_qual = av / reads
               average_size = sizes / reads
            else:
               average_qual = 0
               average_size = 0
        else:
            average_qual = 0
            average_size = 0

        total_reads_counter += reads
        overall_av_counter += average_qual
        overall_size_counter += average_size
        read_list.append(reads)
        av_list.append(average_qual)
        size_list.append(average_size)

    total_reads = total_reads_counter
    total_reads_av = total_reads_counter / len(f)
    overall_av = overall_av_counter / len(f)
    one_direction = total_reads / float(2)
    min_read = min(read_list)
    max_read = max(read_list)
    min_av = min(av_list)
    max_av = max(av_list)
    average_size = overall_size_counter / len(f)
    min_size = min(size_list)
    max_size = max(size_list)

    print("Total reads: %s \
          \nAverage reads per file: %s \
          \nReads into 1 direction: %s \
          \nMin reads: %s \
          \nMax reads: %s \
          \nOverall quality: %s \
          \nWorst average quality: %s \
          \nBest average quality: %s \
          \nAverage read length: %s \
          \nMin length: %s \
          \nMax length: %s" %(total_reads, total_reads_av, one_direction,min_read,max_read,overall_av,min_av,max_av,average_size,min_size,max_size))
