#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
This script sorts all .bam files of a given input directory
"""

import sys
import os
import glob

def aligned_counts(ifile1):
    '''
    count how many alignments are aligned back to genome,
    ifile1 is a sorted bam file
    '''
    import HTSeq
    sortedbamfile= HTSeq.BAM_Reader(ifile1)
    aligned_counts=0
    unaligned_counts=0
    for almnt in sortedbamfile:
        if almnt.aligned:
            aligned_counts+= 1
        else:
            unaligned_counts+=1
    sum = float(aligned_counts) + float(unaligned_counts)
    ratio = (aligned_counts / sum) * float(100)
    #print "number of aligned tags of %s is %d " % (ifile1, aligned_counts)
    #print "number of unaligned tags of %s is %d "% (ifile1, unaligned_counts)
    return ratio

directory = "/home/benflies/NGS_Data/BAM/bam_tst/velona"

bam_file_list = [f for f in glob.iglob(directory+"/*.bam")]

for file in bam_file_list:

    print file

    aligned = aligned_counts(file)

    print "File %s : %s" %(file,aligned)
