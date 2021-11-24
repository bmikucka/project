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
import re

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

   spfeatures = {
      "NULL":      (0),
      'ACT_SITE':  (2**0),
      'BINDING':   (2**1),
      'CA_BIND':   (2**2),
      'DNA_BIND':  (2**3),
      'NP_BIND':   (2**4),
      'METAL':     (2**5),
      'MOD_RES':   (2**6),
      'CARBOHYD':  (2**7),
      'MOTIF':     (2**8),
      'LIPID':     (2**9),
      'DISULFID':  (2**10),
      'CROSSLNK':  (2**11)
   }

   #relevant_fts = ["ACT_SITE", "BINDING", "CA_BIND", "DNA_BIND", "NP_BIND", "METAL", 
   #"MOD_RES", "CARBOHYD", "MOTIF", "LIPID"]
   #special_fts = ["DISULFID", "CROSSLNK"]

   result = 0

   for line in sprot_lines:
      #filter for lines with residue numbers
      feature_line = re.findall("^FT   [A-Z]", line)
      if feature_line:
         #replace multiples of whitespaces
         line = ' '.join(line.split())
         #split the string by white spaces
         info_list = line.split()
         if '..' in info_list[2]:
            #for features with range of residues
            res_range = info_list[2].replace('..', ' ')
            res_range = res_range.split()
            start = int(res_range[0])
            stop = int(res_range[1])
         else:
            #for features at one residue
            if len(info_list) == 3:
               start = int(info_list[2])
               stop = int(info_list[2])
            #for features with residue numbers separated by spaces
            elif len(info_list) == 4:
               start = int(info_list[2])
               stop = int(info_list[3])

         #check for relevant features
         
         if (info_list[1] == 'DISULFID' or info_list[1] == 'CROSSLNK'):
            if (res_of_interest == start and res_of_interest == stop):
               print ('bad')
               print (info_list[1])
               #add the feature to the result (if only mutated features)
               result = result | spfeatures[info_list[1]][0]
            else:
               print ('ok')
            #if nonmutated features also need to be added then put the result part here
         elif info_list[1] in spfeatures:
            #check if residue of interest is in the feature's range
            if (res_of_interest >= start and res_of_interest <= stop):
               #bad for the mutation to be in a relevant feature
               print ('bad')
               print (info_list[1])
               #add the feature to the result (if only mutated features)
               result = result | spfeatures[info_list[1]][0]
            else: 
               print ('ok')
            #if nonmutated features also need to be added then put the result part here
         return result

   #create a string to print binary results
   result = bin(result) 
   #make result into a string 
   result = str(result)
   print (result)

#*************************************************************************






#take file and residue number from command line

sprot_file = sys.argv[1]
res_of_interest = int(sys.argv[2])

check_feature (sprot_file)




   
      
