#!/usr/bin/python3

"""
Program: get_mcc
File:    get_mcc.py

Version:    V1.0
Date:       14.03.22
Function:   List test MCC values from Weka output files


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
file_lines = read_file(sys.argv[1])

result_lines = []

# get the second (so testing not training result) line that starts with Weighted Avg. 
for line in file_lines:
   if line.startswith("Weighted"):
      result_lines.append(line)

result_line = result_lines[1]
result_line = ' '.join(result_line.split())
results = result_line.split()

mcc = results [7]

# print file name and MCC result

print (sys.argv[1], mcc)


# to run
# for file in ./*.out 
# do
# ~/bin/get_mcc.py $file >> mcc_results.txt
# done

