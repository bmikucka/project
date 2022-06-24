#!/bin/python3

#calculate mcc score based on file with TP/TN/FP/FN results across all datapoints for one fold


#*************************************************************************
# Import libraries

import sys
import math
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
results = read_file(sys.argv[1])

TP = 0
FN = 0
FP = 0
TN = 0

for result in results:
   if result == "TP":
      TP = TP + 1
   elif result == "FN":
      FN = FN + 1
   elif result == "FP":
      FP = FP + 1
   elif result == "TN":
      TN = TN + 1


mcc = (TP*TN - FP*FN) / math.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))

print (mcc)