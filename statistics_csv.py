#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
This script takes multiple directories as input,
searches each of them for .fastq.gz files and
writes statistics into a .csv
"""

from __future__ import absolute_import, division, print_function

import sys
import os
import gzip
import csv
import re
import operator
import numpy as n

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

for arg in directories:
    os.chdir("/home")
    path_fastq = arg
    os.chdir(os.path.dirname(os.path.abspath(arg)))
    files = os.listdir(path_fastq)

    output_name = "statistics_2.csv"

    with open(os.path.join(path_fastq,output_name),"w") as f_out:

        #search fastq files:
        f = [os.path.join(root,name)
                     for root, dirs, files in os.walk(path_fastq)
                     for name in files
                     if name.endswith((".fastq", ".fastq.gz"))]

        f.sort()

        collected = {}

        for file in f:

            handle = gzip.open(file)
            reads = 0
            average = []
            av = 0
            sizes = 0

            basename = get_basename(file, ['.fastq', '.fastq.gz'], cut_paired=True)
            name = get_basename(file, ['.fastq', '.fastq.gz'])

            #count reads & calculate average
            for rec in SeqIO.parse(handle, "fastq"):
                reads += 1
                size = len(rec)
                sizes += len(rec)
                if reads % 5000 == 0:
                    print("currently in file %s at read %s" %(name,reads))
                count_qual = 0
                for qual in rec.letter_annotations["phred_quality"]:
                     count_qual += qual
                if len(rec) > 0:
                    average = count_qual / len(rec)
                else:
                    average = 0

                av += average

            if reads > 0:
                if average > 0:
                    average_qual = av / reads
                    av_size = sizes / reads
                else:
                    average_qual = 0
                    av_size = 0
            else:
                average_qual = 0
                av_size = 0

            # collect data
            data = [file, reads, average_qual,av_size]

            if basename not in collected:
                collected[basename] = {
                    "reads": reads,
                    "quals": [average_qual],
                    "length": [av_size]
                }
            else:
                 collected[basename]['quals'].append(average_qual)
                 collected[basename]['length'].append(av_size)
                 collected[basename]['average qual'] = float(sum(collected[basename]['quals'])) / float(2)
                 collected[basename]['average size'] = float(sum(collected[basename]['length'])) / float(2)

        # ouput data into a csv
        writer = csv.writer(f_out)
        writer.writerow(['sample', 'reads', 'forward quality', 'reverse quality', 'overall average quality', 'forward length', 'reverse length', 'average length' ])
        for sample, values in collected.items():
            data = [sample, values['reads'], values['quals'][0], values['quals'][1], values['average qual'], values['length'][0], values['length'][1], values['average size']]
	    writer.writerow(data)
