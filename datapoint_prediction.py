#!/bin/python3

#get TN/TP/FN/FP result from confidence scores

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

first_line_elements = file_lines[0].split()
actual = first_line_elements[0]
scores = []

for line in file_lines:

   answers = line.split()
   score = answers[1]
   #make list of confidence scores across all balancing runs to get average 
   scores.append(score)

#have list of scores to average and the real classification

avg_score = sum(scores) / len(scores)
if avg_score > 0.5:
   predicted = "1:PD"
elif avg_score < 0.5:
   predicted = "2:SNP"


if actual == "1:PD" and predicted == "1:PD":
   result = "TP"
elif actual == "1:PD" and predicted == "2:SNP":
   result = "FN"
elif actual == "2:SNP" and predicted == "1:PD":
   result = "FP"
elif actual == "2:SNP" and predicted == "2:SNP":
   result = "TN"


print (result)