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
def format_line_range (file_line):
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
   min_res = int(info_list [2])
   max_res = int(info_list [3])

   return (ID, feature, min_res, max_res)

#*************************************************************************
def format_line_single (file_line):
   """ Formats a line of the SwissProt file to split into ID, feature and 
   one residue that constitutes the feature.

   Input:   file_line   --- line of SwissProt file
   Return:  ID          --- line ID (FT if a feature)
            feature     --- type of feature
            res         --- residue number at which the feature starts
            
   16.11.21    Original    By: BAM

   """
    
   line = ' '.join(file_line.split())
   #split the string by white spaces
   info_list = line.split()
   ID = info_list [0]
   feature = info_list [1]
   res = int(info_list [2])

   return (ID, feature, res)

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

   for line in sprot_lines:
      #filter out lines with FT and no / - this doesnt work
      if (line.startswith('FT')) and ('/' not in line):
         #for features across a range or residues
         if '..' in line:
            (ID, feature, min_res, max_res) = format_line_range (line)
            #check if feature is relevant
            if feature in relevant_fts:
               #check if residue of interest it in the range of feature
               if (residue_of_interest >= min_res and residue_of_interest <= max_res):
                  #bad for the mutation to be in a relevant feature
                  print ('bad')
                  print (feature)
               else: print ('ok')

         #for features at one residue
         else:
            (ID, feature, res) = format_line_single (line)
            #chech if feature is relevant
            if feature in relevant_fts:
               #check if residue of interest is the residue of the feature
               if residue_of_interest == res:
                  #bad for the mutation to be in a relevant feature
                  print ("bad")
                  print (feature)
               else: print('ok')
               #need to change this so it doesnt print ok/bad separately 
               #for the single residue features


#*************************************************************************


#take file and residue number from command line

sprot_file = sys.argv[1]
residue_of_interest = int(sys.argv[2])

check_feature (sprot_file)


