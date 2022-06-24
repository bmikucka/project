#!/bin/python3

#get TN/TP/FN/FP result from the file and print

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
   if line.startswith("     a"):
      title = file_lines.index(line)
      title = int(title)
      #get index of lines of confusion matrix
      positive_line_index = title + 1
      negative_line_index = title + 2

      #get lines of confusion matrix
      positive_line = file_lines[positive_line_index]
      negative_line = file_lines[negative_line_index]

      #make list of numbers classified as positives and negatives
      positives = positive_line.split()
      negatives = negative_line.split()

      #get True and False Values

      TP = positives[0]
      FN = positives[1]
      FP = negatives[0]
      TN = negatives[1]

      #check which on it is

      if TP == "1":
         result = "TP"
      elif FN == "1":
         result = "FN"
      elif FP == "1":
         result = "FP"
      elif TN == "1":
         result = "TN"


print (result)