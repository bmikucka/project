#!/bin/python3

#get confidence score for one datapoint in one balancing run

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


for line in file_lines:
   if line.startswith("        1"):

      answers = line.split()

      score = answers[-1]
      score = float(prediction)
      score = score/2

      predicted = answers[2]

      actual = answers[1]

      if predicted == "1:PD":
         result = 0.5 + score
      elif predicted == "2:SNP":
         result = 0.5 - score

print (actual, result)