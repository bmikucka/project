#!/usr/bin/python3

"""
Program: mcc_avg
File:    mcc_avg.py

Version:    V1.0
Date:       14.03.2022
Function:   Calculates average and stdev values for mcc_results files.


Author: Barbara A. Mikucka

--------------------------------------------------------------------------
Usage:
======

--------------------------------------------------------------------------
Revision history:
=================
V1.0  14.03.2022    Original    By: BAM
"""

#*************************************************************************
# Import libraries

import sys
import statistics

#*************************************************************************
def read_file (filename):
   """ Read input file and split by lines into a list.

   Input:   filename    --- Swiss Prot File 
   Return:  file_lines  --- list of lines from the SwissProt file 

   16.11.21    Original    By: BAM

   """
   with open (filename) as file:
      #read all lines
      file_lines = file.readlines()
      return file_lines

#*************************************************************************

# read file as lines
file_lines = read_file(sys.argv[1])

# compile MCC results in list
mcc_results = []

for line in file_lines:
   line_infos = line.split()
   mcc = float(line_infos[1])
   mcc_results.append(mcc)

# calculate average
mcc_avg = (sum(mcc_results)) / (len(mcc_results))
# calculate standard deviation
st_dev = statistics.stdev(mcc_results)

print (mcc_avg, st_dev)
