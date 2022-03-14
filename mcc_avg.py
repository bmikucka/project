#!/usr/bin/python3

import sys

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


file_lines = read_file(sys.argv[1])

mcc_results = []
for line in file_lines:
   line_infos = line.split()
   mcc_int = int(line_infos[1])
   mcc_results.append(mcc_int)

mcc_avg = int(sum(mcc_results)) / int(len(mcc_results))

print (mcc_avg)
