#!/usr/bin/python3

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
   lines = read_file (filename)
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
   lines = read_file (filename_a)
   for line in lines:
      infos = line.split()
      mcc_results_a.append(infos[1])

   mcc_results_b = []    #list of MCC results from dataset B
   lines = read_file (filename_b)
   for line in lines:
      infos = line.split()
      mcc_results_a.append(infos[1])



   # Getting data in to a dictionary
   data = {'A': mcc_results_a,
            'B': mcc_results_b}

   # Dictionary to Dataframe
   df = pd.DataFrame(data)

   #check if difference is significant
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




file_a = sys.argv[1]
file_b = sys.argv[2]

print("\n \nDATASET A:")
normal_test(file_a)

print("\nDATASET B")
normal_test(file_b)

print("\n \nMann-Whitney U Test:")
mwu_test(file_a, file_b)