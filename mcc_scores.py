#!/bin/python3

#get mean MCC value for the ten folds results

#*************************************************************************
# Import libraries

import sys
import statistics
#*************************************************************************
def read_file (filename):
   """ Read input file and split by lines into a list.

   Input:   filename    
   Return:  file_lines  --- list of lines from the input file 

   16.11.21    Original    By: BAM

   """
   with open (filename) as file:
      #read all lines
      file_lines = file.readlines()
      return file_lines

#*************************************************************************

# read file
lines = read_file(sys.argv[1])

mccs = []

for line in lines:
   if line != "":
      mccs.append(line)

mcc_avg = statistics.mean(mccs)

print (mcc_avg)


