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
import os
import json
from urllib.request import urlopen
import atomium

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
def process_resnum (res_id):
   """Separate chain id, residue number and insert from the resnum input

   Input:  res_id     --- input from command line: [chain]resnum[insert]
   Return: chain_id   --- Optional, chain id of the PDB file 
           resnum_pdb --- PDB number of the mutant residue 
           insert     --- insert code


   14.01.22    Original    By: BAM

   """

   #if the first character is a letter - make that the chain_id
   if (res_id[0]).isalpha():
      chain_id = res_id[0]
      res_id = res_id [1:]
   else: chain_id = ''

   if (res_id[-1]).isalpha():
      insert = res_id [-1]
      res_id = res_id [:-1]
   else: insert = ''

   resnum_pdb = res_id

   return (chain_id, resnum_pdb, insert)


#*************************************************************************
def get_pdb_code (pdbfile):
   """Get PDB code from PDB file name

   Input:  pdbfile    --- PDB file
   Return: pdb_code   --- PDB code 

   14.01.22    Original    By: BAM

   """
   
   #file name as string
   pdb_name = str(pdbfile)
   #take out the pdb code
   pdb_code = pdb_name[18:-4]

   return pdb_code


#*************************************************************************
def read_url_sprot (sprot_ac):

   """Read SwissProt file from the web.
   
   Input:  sprot_ac   --- SwissProt accession number
   Return: read file  
   

   14.01.22    Original    By: BAM

   """

   url = 'https://www.uniprot.org/uniprot/{}.txt'.format(sprot_ac)
   file = urlopen(url)
   return (file.read())


#*************************************************************************
def get_sprot_str (sprot_ac):
   """Read SwissProt file from the web.

   Input:  sprot_ac   --- SwissProt accession number
   Return: sprot_str  --- SwissProt file as string 

   """

   sprot_file_byt = read_url_sprot (sprot_ac)
   #this is now a byte object - convert into string
   encoding = 'utf-8'
   sprot_str = sprot_file_byt.decode(encoding)

   return sprot_str


#*************************************************************************
def get_ft_residues (sprot_str, res_of_interest):
   """ Returns list of residues that are in relevant features based on 
   SwissProt file.

   Input:   sprot_file      --- Swiss Prot File   
            res_of_interest --- Number of residue being checked
   Output:  sp_ft_residues  --- List of lists (one per feature) of the 
                                SwissProt residue numbers in that feature
           

   09.12.21    Original    By: BAM

   """
   #get list of lines from the file
 
   #make list of relevant features


   #list of residues in relevant features
   sp_ft_residues = []

   #dictionary of relevant features
   spfeatures = ['ACT_SITE', 
      'BINDING', 
      'CA_BIND', 
      'DNA_BIND', 
      'NP_BIND', 
      'METAL', 
      'MOD_RES', 
      'CARBOHYD', 
      'MOTIF', 
      'LIPID']
   


   #split sprot_file into lines
   sprot_lines = sprot_str.split('\n')

   for line in sprot_lines:
      #filter for lines with residue numbers
      feature_line = re.findall("^FT   [A-Z]", line)

      if feature_line:
         #replace multiples of whitespaces
         line = ' '.join(line.split())
         #split the string by white spaces
         info_list = line.split()

         if info_list[1] in spfeatures:
            if '..' in info_list [2]:
               #for features with a range of residues
               res_range_str = info_list[2].replace('..', ' ')
               two_residues = res_range_str.split() 
               start = int(two_residues[0])
               stop = int(two_residues[1])
               #get numbers of all the residues into a list
               residues = range(start, stop)
               #add everything from start to stop to the sp_ft_residues list
               sp_ft_residues.append(residues)

            else:
               #for one residue long features
               start = int(info_list[2])
               #make into a list
               residues = [start]
               #add the residue to the list
               sp_ft_residues.append(residues)
    
   return sp_ft_residues 


#*************************************************************************
def pdb_ft_list (sp_ft_residues, sprot_ac, pdb_code, chain_id):
   """ Translates the SwissProt numbering of residues (in list) to PDB.

   Input:   sp_ft_residues  --- List of lists (one per feature) of 
                                SwissProt residue numbers.
   Output:  pdb_ft_residues --- List in same format as input but with
                                PDB numbering.


   14.01.22    Original    By: BAM

   """

   pdb_ft_residues = []
   #list of lists with residue numbers in PDB numbering for each feature
   for feature in sp_ft_residues:
      feature_residues_pdb = []
      for residue in feature:
         residue_pdb = sws_pdb(sprot_ac, residue, pdb_code, chain_id)
         feature_residues_pdb.append(residue_pdb)
      pdb_ft_residues.append(feature_residues_pdb)

   #remove empty strings
   temp_list_fts = []
   for feature in pdb_ft_residues:
      temp_list_res = []
      for number in feature:
         if number != '' and number != None: 
            temp_list_res.append(number)
            #list for each feature including non empty residue numbers
      numbers = temp_list_res

      #only add back nonempty lists
      if numbers != []:
         temp_list_fts.append(numbers)

   pdb_ft_residues = temp_list_fts

   return pdb_ft_residues


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
def read_url_swspdb (uniprot_ac, uniprot_resid):
   """Reads text file from url PDBSWS to get information on PDB files 
   related to the UniProt accession number.

   Input:   uniprot_ac      --- UniProt accession number  
            uniprot_resid   --- Number for residue of interest in UniProt file
   Output:  file            --- Text file from url as a byt
   
   09.12.21    Original    By: BAM

   """
  
   url = 'http://www.bioinf.org.uk/servers/pdbsws/query.cgi?plain=1&qtype=ac&id={}&res={}'.format(uniprot_ac, uniprot_resid)
   file = urlopen(url)
   return (file.read())


#*************************************************************************
def sws_pdb (uniprot_ac, uniprot_resid, pdb_code, chain_id):
   """ Returns PDB format residue number for the UniProt residue number.

   Input:   uniprot_ac      --- UniProt/SwissProt accession number  
            uniprot_resid   --- Number for residue of interest in UniProt/ 
                                SwissProt file (integer)
            pdb_code        --- The relevant number will be returned for 
                                this PDB code
            chain_id        --- The relevant number will be returned for 
                                this chain_id
   Output:  resnum_pdb      --- number of the residue in the PDB file
   
   09.12.21    Original    By: BAM

   """

   #empty for all SwissProt residues that do not have a PDB equivalent
   resnum_pdb = ''

   #change integer into string
   str(uniprot_resid)

   pdb_file_byt = read_url_swspdb(uniprot_ac, uniprot_resid)
   #this is now a byte object - convert into string
   encoding = 'utf-8'
   pdb_file_str = pdb_file_byt.decode(encoding)

   #split into list of parts for each PDB chain
   pdb_file_parts = pdb_file_str.split('//')

   for file_part in pdb_file_parts:
      if file_part.startswith('\n'):
         file_part = file_part[1:]
      lines = file_part.split('\n')
      #remove empty lines

      if lines != ['']:
         nonempty_lines = lines

         pdb_line = nonempty_lines [0]
         chain_line = nonempty_lines [1]

         if (pdb_line [5:9] == pdb_code) and (chain_line [7] == chain_id):
            residue_line = nonempty_lines [2]
            resnum_pdb = residue_line [7:]
            return resnum_pdb


#*************************************************************************
def feature_distance (pdb_code, chain_id, resnum_pdb, pdb_ft_residues):
   """ Processes lists of residue numbers, gets shortest distance 
   between a residue of interest atom and one of the feature atoms.

   Input:   pdb_code          --- PDB code 
            chain_id          --- Chain ID
            resnum_pdb        --- residue of interest (mutant) number 
            pdb_ft_residues   --- list of lists (for each feature) with 
                                   residue numbers in that feature
   Output:  feature_distances --- List of lists (one for each feature)
                                  with distance to the mutant residue 
                                  and the closest feature residue number

   
   10.12.21    Original    By: BAM

   """

   #list of info for each feature as output:
   feature_distances= []

   #get list of atoms in the mutant residue
   model = atomium.fetch(pdb_code).model
   chain = model.chain(chain_id)
   residue = chain.residue(f"{chain_id}.{resnum_pdb}")

   atoms = []
   for atom in residue.atoms():
      atoms.append(atom.id)



   #unpack the list of features 

   #for each of the features
   for residue_list in pdb_ft_residues:
      min_dist = 1000 #reset the smallest distance after each feature

      #for each of the residues in that feature
      for aa in residue_list:
         #get list of atoms in that residue
         ft_atoms = []
         residue = chain.residue(f"{chain_id}.{aa}")
         for atom in residue.atoms():
            ft_atoms.append(atom.id) 

         #get distance between each atom combination
         for x in ft_atoms:
            for y in atoms:
               d = model.atom(int(x)).distance_to(model.atom(int(y)))
               #record the smallest distance for each feature and the  
               #corresponding residue number from the feature
               if d < min_dist:
                  min_dist = d
                  closest_res = aa

      #for each feature add a list with the recorded information
      feature_distances.append ([min_dist, closest_res]) 

   return (feature_distances)


#*************************************************************************
def read_url_pdbsws (pdb, chain, pdb_residue):
   """ Get UniProt id and residue number equivalents from a specified PDB 
   residue.

   Input:   pdb            --- 4 letter PDB code
            chain          --- Chain ID
            pdb_residue    --- PDB Residue number
    Output:  file          --- Text file from url as a byt

   11.12.21    Original    By: BAM

   """
   url = 'http://www.bioinf.org.uk/servers/pdbsws/query.cgi?plain=1&qtype=pdb&id={}&chain={}&res={}'.format(pdb, chain, pdb_residue)
   file = urlopen(url)
   return (file.read())


#*************************************************************************
def pdb_sws (pdb, chain, pdb_residue):
   """ Get UniProt accession number and residue number equivalents from a 
   specified PDB residue.

   Input:   pdb            --- 4 letter PDB code
            chain          --- Chain ID
            pdb_residue    --- PDB Residue number
   Output:  sprot_residue  --- SwissProt residue number
            sprot_ac       --- SwissProt AC

   11.12.21    Original    By: BAM

   """

   info_byt = read_url_pdbsws (pdb, chain, pdb_residue)
   #this is now a byte object - convert into string
   encoding = 'utf-8'
   info_str = info_byt.decode(encoding)

   #split into lines
   info_lines = info_str.split('\n')

   #get the line that starts with AC: and UPCOUNT:
   for line in info_lines:
      if line.startswith('AC:'):
         words = line.split(' ')
         #get the second word in the line - the accession number
         sprot_ac = words[1]
      if line.startswith('UPCOUNT:'):
         words = line.split(' ')
         #get the residue number
         sprot_residue = words [1]

   return (sprot_ac, sprot_residue)


#*************************************************************************
def residue_to_feature (res_id, sprot_ac):
  """ Gets the ID of feature that the provided residue number is in.

   Input:   res_id   --- residue number from SwissProt file
            sprot_ac --- SwissProt accession number
   Output:  feature  --- SwissProt feature ID

   13.12.21    Original    By: BAM

   """ 
  

  #dictionary of relevant features
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



  #get the SwissProt file
  file = read_url_sprot (sprot_ac)

  encoding = 'utf-8'
  #make into a string
  sprot_str = file.decode(encoding)

  #split text into lines
  sprot_lines = sprot_str.split('\n')


  for line in sprot_lines:
   #filter for lines with features
   feature_line = re.findall("^FT   [A-Z]", line)

   if feature_line:
      #replace multiples of whitespaces
      line = ' '.join(line.split())
      #split the string by white spaces
      info_list = line.split()
      
      #check the lines with features of interest
      if info_list[1] in spfeatures:

         #if the residues were separated by a space
         if len(info_list) == 4:
            #if a 2 residue feature
            if (info_list[1] == 'DISULFID' or info_list[1] == 'CROSSLNK'):
               if res_id == info_list[2] or res_id == info_list[3]:
                  feature = info_list[1]
                  return(feature)

            #if a range residue feature
            else: 
               residue_numbers = range(int(info_list[2]), int(info_list[3]))
               if int(res_id) in residue_rumbers:
                  feature = info_list[1]
                  return(feature)

         #if it's a 1 residue feature or were split by .. and not a space
         elif len(info_list) == 3:

            #for multiple residue features
            if '..' in info_list[2]:
               
               number_str = info_list[2].replace('..', ' ')
               #number is a string of 2 number separated by a space
               residues = number_str.split()
               #residues is a list of the two residues in number
               #now have a list of two residue numbers

               #if a 2 residue feature
               if (info_list[1] == 'DISULFID' or info_list[1] == 'CROSSLNK'):

                  if res_id in residues:
                     feature = info_list[1]
                     return (feature)

               #if a feature across the range of residues
               else:
                  residue_numbers = range(int(residues[0]), int(resdiues[1]))
                  if int(res_id) in residue_rumbers:
                     feature = info_list[1]
                     return(feature)

            #for single residue features
            else:
               if info_list[2] == res_id:
                  feature = info_list[1]
                  return(feature)

         
#*************************************************************************
def get_bool_results (distances_list):
   """ Get result BAD or OK.

   Input:   final_infos  --- dictionary 
   Output:  result       --- BAD - likely to have an effect or OK - 
                             unlikely to have an effect

   14.01.22    Original    By: BAM

   """

   for d in distances_list:
      if d < 4:
         result = "BAD"
         return result
      else: result = "OK"
   
   return result


#*************************************************************************
def check_cache (res_id, newaa, pdbfile, opts):
   """ Check if results have been cached.

   Input:   res_id        --- input residue 
            newaa         --- input new mutant residue
            pdbfile       --- input PDB file
            opts          --- list of optional command line arguments
   Output:  cached_file   --- result of the program in json format, if 
                              cached before, empty string if not cached

   18.01.22    Original    By: BAM

   """


   #make file name using the input
   file_name = ("{}_{}_{}.txt").format(pdbfile, res_id, newaa)
   file_name = file_name.replace("/", "_")
   path = './cachedir/SprotFTdist'

   full_name = os.path.join(path, file_name)
   
   #if -f in opts return empty string (force to run)
   if '-f' in opts or '-force' in opts:
      return ('')

   #check if filename in directory

   #this doesnt work
   if os.path.isfile(full_name):
      #if in the directory then return the content
      with open (full_name, 'r') as file:
         cache_file = file.read()
         return cache_file

   #if not in the directory then restun empty string
   else: 
      return ('')

#*************************************************************************
def write_cache (res_id, newaa, pdbfile, output_str):
   """ Write cache file for this input.

   Input:   res_id   --- input residue 
            newaa    --- input new mutant residue
            pdbfile  --- input PDB file
            output   --- content of the new file

   18.01.22    Original    By: BAM

   """
   file_name = ("{}_{}_{}.txt").format(pdbfile, res_id, newaa)
   file_name = file_name.replace("/", "_")
   path = './cachedir/SprotFTdist'

   #convert output string into byte
   output_encoded = bytes(output_str,'UTF-8')

   with open(os.path.join(path, file_name), 'wb') as cache_file:
    cache_file.write(output_encoded)

   return


#*************************************************************************
