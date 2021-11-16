"""
Program: sprotfeatures
File:    sprotfeatures.py

Version:    V1.0
Date:       16.11.2021
Function:   Identifies if residue of interest is marked with a relevant 
            feature in a SwissProt file


Author: Barbara A. Mikucka

--------------------------------------------------------------------------
Description:


--------------------------------------------------------------------------
Usage:
======

--------------------------------------------------------------------------
Revision history:
=================
V1.0  16.11.21    Original    By: BAM
"""

#*************************************************************************
# Import libraries

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
def format_line (file_line):
   """ Formats a line of the SwissProt file to split into ID, feature and 
   range of residues that constitute the feature.

   Input:   file_line   --- line of SwissProt file
   Return:  ID          --- line ID (FT if a feature)
            feature     --- type of feature
            min_res     --- residue number at which the feature starts
            max_res     --- residue number at which the feature ends

   16.11.21    Original    By: BAM

   """
   #change .. for a white space
   line = file_line.replace('..', ' ')
    
   line = ' '.join(line.split())
   #split the string by white spaces
   info_list = line.split()
   ID = info_list [0]
   feature = info_list [1]
   min_res = info_list [2]
   #problem - doesnt always have two residues
   max_res = info_list [3]

   return (ID, feature, min_res, max_res)

#*************************************************************************
def check_feature (sprot_file):
   """ Checks if the line refers to a feature, if the feature is relevant 
   and if the residue of interest is in the range of the feature. Prints 
   "ok" if the residue is not in a feature and "bad" if within a feature.

   Input:   sprot_file  --- Swiss Prot File      

   16.11.21    Original    By: BAM

   """
   #get list of lines from the file
   sprot_lines = read_file(sprot_file)
   #make list of relevant features
   relevant_fts = ["ACT_SITE", "BINDING", "CA_BIND", "DNA_BIND", "NP_BIND", "METAL", 
   "MOD_RES", "CARBOHYD", "MOTIF", "LIPID", "DISULFID", "CROSSLNK"]

   #check if line with a feature
   for line in sprot_lines:
      if line.startswith('FT'):
         (ID, feature, min_res, max_res) = format_line (line)
         #PROBLEM - not all lines have 4 strings in them
         #check if the feature is relevant
         if feature in relevant_fts:
            #check if the residue of interest is in the residue range
            if (residue_of_interest >= min_res and residue_of_interest <= max_res):
               #bad for mutation in residue of interest if within a relevant feature
               print ("bad")
               print (feature)
            else: print ('ok')




#*************************************************************************


#take file and residue number from command line

sprot_file = sys.argv[1]
residue_of_interest = sys.argv[2]

check_feature (sprot_file)


