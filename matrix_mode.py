#!/bin/python3

#take file of TP/FP/TN/FN and get the mode and print

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
results = read_file(sys.argv[1])

mode = statistics.mode(results)

print (mode)