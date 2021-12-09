"""
Program: sprotfeatures_functions
File:    sprotfeatures_functions.py

Version:    V1.0
Date:       09.12.2021
Function:   File of funcitons for sprotfeatures


Author: Barbara A. Mikucka

--------------------------------------------------------------------------
Description:


--------------------------------------------------------------------------
Usage:
======

--------------------------------------------------------------------------
Revision history:
=================
V1.0  09.12.21    Original    By: BAM
"""

#*************************************************************************
# Import libraries

import sys
import re
from urllib.request import urlopen

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
def check_feature (sprot_file, res_of_interest):
   """ Checks if the line refers to a feature, if the feature is relevant 
   and if the residue of interest is in the range of the feature. Prints 
   "ok" if the residue is not in a feature and "bad" if within a feature.

   Input:   sprot_file      --- Swiss Prot File   
            res_of_interest --- Number of residue being checked
   Output:  outcome         --- 'bad' if residue of interest is in a 
                                feature of interest, otherwise 'ok'
            result          --- string with binary representation of 
                                relevant features the residue is in
            mut_features    --- list of features affected if residue of 
                                interest is mutated
           

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

   result = 0
   #variable to check if the residue of interest is in any relevant features
   a = 0
   #list of features affected by the mutation at the residue
   mut_features = []


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
            res_range_str = info_list[2].replace('..', ' ')
            res_range = res_range_str.split() 
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
            result = result | spfeatures[info_list[1]]

            if (res_of_interest == start or res_of_interest == stop):
               a = a + 1
               #add the feature to the list of affected features
               mut_features.append(info_list[1])

         elif info_list[1] in spfeatures:
            result = result | spfeatures[info_list[1]]

            #check if residue of interest is in the feature's range
            if (res_of_interest >= start and res_of_interest <= stop):
               #bad for the mutation to be in a relevant feature
               a = a + 1
               #add the feature to the list of affected features
               mut_features.append(info_list[1])
            


   #create a string to print binary results
   result = bin(result)[2:].zfill(12)
   #make result into a string 
   result = str(result)

   if a == 0:
      outcome = 'ok'
   else: 
      outcome = 'bad'

   return (outcome, mut_features, result)



#*************************************************************************
def read_url_file (uniprot_ac, uniprot_resid):
   """Reads text file from url PDBSWS to get information on PDB files 
   related to the UniProt accession number.

   Input:   uniprot_ac      --- UniProt accession number  
            uniprot_resid   --- Number for residue of interest in UniProt file
   Output:  file            --- Text file from url as a string
   
   09.12.21    Original    By: BAM

   """
  
   #get PDB codes corresponding to the UniProt accession number using PDBSWS
   url = 'http://www.bioinf.org.uk/servers/pdbsws/query.cgi?plain=1&qtype=ac&id={}&res={}'.format(uniprot_ac, uniprot_resid)
   file = urlopen(url)
   return (file.read())


#*************************************************************************

def pdb_sws (uniprot_ac, uniprot_resid):
   """ Returns PDB code, chain and residue number for all PBDs relevant to 
   the given UniProt accession number.

   Input:   uniprot_ac      --- UniProt accession number  
            uniprot_resid   --- Number for residue of interest in UniProt file
   Output:  pdb_infos       --- list of PDB code, chain and resid relevant   
                                to the UniProt accession number 
   
   09.12.21    Original    By: BAM

   """

   #make list of PDB codes
   pdb_infos = []
   pdb_code = ''
   chain = ''
   resid = ''

   pdb_file_byt = read_url_file(uniprot_ac, uniprot_resid)
   #this is now a byte object - convert into string
   encoding = 'utf-8'
   pdb_file_str = pdb_file_byt.decode(encoding)

   #split into list of parts for each PDB chain
   pdb_file_parts = pdb_file_str.split('//')

   #for each part of the file read each line
   for file_part in pdb_file_parts:
      #split into lines
      line = file_part.split('\n')
      
      #from each PDB info part get the pdb code, chain, and resid
      for x in line:
         if x.startswith('PDB:'):
            #split the line into individual words
            words = x.split(' ')
            #thefirst word is PDB:, second is the actual code
            pdb_code = words[1]
            #repeat for chain and resid
         elif x.startswith('CHAIN:'):
            words = x.split(' ')
            chain = words[1]
         elif x.startswith('RESID'):
            words = x.split(' ')
            resid = words[1]

      #get the 3 strings into a list, these represent info on one PDB chain
      infos = [pdb_code, chain, resid]

      #add the list to the list of infos
      pdb_infos.append(infos)

   return pdb_infos


#*************************************************************************

