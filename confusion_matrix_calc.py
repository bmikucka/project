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
   if line.startswith("        1"):

      answers = line.split()

      actual = answers[1]
      predicted = answers[2]

      if actual == "1:PD" and predicted == "1:PD":
         result = "TP"
      elif actual == "1:PD" and predicted == "2:SNP":
         result = "FN"
      elif actual == "2:SNP" and predicted == "1:PD":
         result = "FP"
      elif actual == "2:SNP" and predicted == "2:SNP":
         result = "TN"


print (result)