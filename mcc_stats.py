#!/usr/bin/python3

"""
Program: mcc_stats
File:    mcc_stats.py

Version:    V1.0
Date:       16.04.2022
Function:   Does Shaphiro-Wilk and Mann-Whitney U tests on two datasets.


Author: Barbara A. Mikucka

--------------------------------------------------------------------------
Usage:
======

--------------------------------------------------------------------------
Revision history:
=================
V1.0  16.04.2022    Original    By: BAM
"""

#*************************************************************************
# Import libraries

from scipy import stats
from scipy.stats import mannwhitneyu
import sys
import pandas as pd


#*************************************************************************
def read_file (filename):
   """ Read input file and split by lines into a list.

   Input:   filename    --- Swiss Prot File 
   Return:  file_lines  --- list of lines from the file 

   15.04.22    Original    By: BAM

   """
   with open (filename) as file:
      #read all lines
      file_lines = file.readlines()
      return file_lines

#*************************************************************************
def normal_test (filename):
   """ Test the dataset is normally distributed.

   Input:   filename  --- list of datatpoints

   15.04.22    Original    By: BAM

   """

   mcc_results = []    #list of MCC results from dataset 
   # read file as lines
   lines = read_file (filename)
   # get mcc values from file lines
   for line in lines:
      infos = line.split()
      mcc_results.append(infos[1])

   #check if data in dataset A follows normal distribution
   shapiro_test = stats.shapiro(mcc_results)

   print("p statistic is:")
   print (shapiro_test.pvalue)
   print("W statistic is:")
   print(shapiro_test.statistic)

   if shapiro_test.pvalue < 0.05:
      print("The data is NOT normally distributed")

   else:
      print ("The data is normally distributed")

   return 

#*************************************************************************
def mwu_test (filename_a, filename_b):
   """ Performs Mann-Whitney U Test on two datatsets.

   Input:   filename_a  --- file with dataset a
            filename_b  --- file with dataset b

   15.04.22    Original    By: BAM

   """

   mcc_results_a = []    #list of MCC results from dataset A
   # read file as lines
   lines = read_file (filename_a)
   # get mcc values from lines for dataset A
   for line in lines:
      infos = line.split()
      mcc_results_a.append(infos[1])

   mcc_results_b = []    #list of MCC results from dataset B
   # read file as lines
   lines = read_file (filename_b)
   #get mcc values from lines for dataset B
   for line in lines:
      infos = line.split()
      mcc_results_b.append(infos[1])



   # Getting data in a dictionary format
   data = {'A': mcc_results_a,
            'B': mcc_results_b}

   # Dictionary to Dataframe
   df = pd.DataFrame(data)

   #check if difference is significant - MWU TEST
   U1, p = stats.mannwhitneyu(df['A'], df['B'])

   print("p value is:")
   print (p)
   print("U value is:")
   print(U1)

   if p < 0.05:
      print("The difference is statistically significant")

   else:
      print ("The difference is NOT statistically significant")

   return 

#*************************************************************************

# files from command line
file_a = sys.argv[1]
file_b = sys.argv[2]

print("\n \nDATASET A:")
# check if dataset A follows normal distribution
normal_test(file_a)

print("\nDATASET B")
# check if dataset B follows normal distribution
normal_test(file_b)

print("\n \nMann-Whitney U Test:")
# run Mann-Whitney U Test
mwu_test(file_a, file_b)