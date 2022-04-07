#!/usr/bin/python3

# get MCC result

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

# get the second line that starts with Weighted Avg. 
for line in file_lines:
   if line.startswith("Weighted"):
      result_lines.append(line)

result_line = result_lines[1]
result_line = ' '.join(result_line.split())
results = result_line.split()

mcc = results [7]

# print file name and MCC result

print  (mcc)


# to run
# for file in ./*.out 
# do
# ~/bin/get_mcc.py $file >> mcc_results.txt
# done

