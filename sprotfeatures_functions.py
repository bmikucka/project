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
def read_url_sprot (sprot_id):

   url = 'https://www.uniprot.org/uniprot/{}.txt'.format(sprot_id)
   file = urlopen(url)
   return (file.read())


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
            if info_list[1] == 'DISULFID' or 'CROSSLNK':
               if '..' in info_list[2]:
                  #for features with range of residues
                  res_range_str = info_list[2].replace('..', ' ')
                  res_range = res_range_str.split() 
                  start = int(res_range[0])
                  stop = int(res_range[1])
                  residues = [start, stop]
                  #only the specific residues are involved
                  sp_ft_residues.append(residues)

                  
               elif len(info_list) == 4:
                  #for features with residue numbers separated by spaces
                  start = int(info_list[2])
                  stop = int(info_list[3])
                  #add the two residues to the sp_ft_residues list in their own list
                  residues = [start, stop]
                  #only the specific residues are involved
                  sp_ft_residues.append(residues)


            else:
               #for features with range of residues
               if '..' in info_list[2]:
                  res_range_str = info_list[2].replace('..', ' ')
                  res_range = res_range_str.split() 
                  start = int(res_range[0])
                  stop = int(res_range[1])
                  #get numbers of all the residues into a list
                  residues = range(start, stop)
                  #add everything from start to stop to the sp_ft_residues list
                  sp_ft_residues.append(residues)
                  
               else:

                  #for features at one residue
                  if len(info_list) == 3:
                     start = int(info_list[2])
                     #add the residue to the sp_ft_residues list
                     sp_ft_residues.append(start)

                  #for features with residue numbers separated by spaces
                  elif len(info_list) == 4:
                     start = int(info_list[2])
                     stop = int(info_list[3])
                     #get numbers of all the residues
                     residues = range(start, stop)
                     #add everything from start to stop to the sp_ft_residues list
                     sp_ft_residues.append(residues)
    
   return sp_ft_residues 




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

def sws_pdb (uniprot_ac, uniprot_resid):
   """ Returns PDB code, chain and residue number for all PBDs relevant to 
   the given UniProt accession number.

   Input:   uniprot_ac      --- UniProt accession number  
            uniprot_resid   --- Number for residue of interest in UniProt file (integer)
   Output:  pdb_infos       --- list of PDB code, chain and resid relevant   
                                to the UniProt accession number 
   
   09.12.21    Original    By: BAM

   """

   #make list of PDB codes
   pdb_infos = []
   pdb_code = ''
   chain = ''
   resid = ''

   #change integer into string
   str(uniprot_resid)

   pdb_file_byt = read_url_swspdb(uniprot_ac, uniprot_resid)
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

      #add the list to the list of infos but only if not empty
      if infos != ['', '', '']:
         pdb_infos.append(infos)

   return pdb_infos


#*************************************************************************
def feature_distance (pdb_infos_res, pdb_infos_fts):
   """ Processes lists of information from PDB about residue of interest 
   and features, gets shortest distance between a residue of interest atom
   and one of the feature atoms.

   Input:   pdb_infos_res     --- list of lists with 3 items (PDB code, chain
                                  identifier, and residue number) for the residue 
                                  of interest
            pdb_infos_fts     --- list of lists for each reasidue with lists with 
                                  3 items (PDB code, chain identifier, and residue 
                                  number) for the resdiues in a feature
   Output:  feature_info_pdb  --- List with an item for each feature from the 
                                  pdb_infos_fts list. Each item is a list 
                                  with the smallest distance between the feature
                                  and the residue of interest, the pdb and chain 
                                  identifiers, and the residue number

   
   10.12.21    Original    By: BAM

   """

   #list of info for each feature as output:
   feature_infos_pdb = []

   #unpack the list of features
   for feature_info in pdb_infos_fts:
      for residue_info in feature_info:
         min_dist = 1000 #reset the smallest distance after each feature
         #for every feature in the protein
         for x in residue_info:
            #for every residue in the feature
            pdb_code = x[0]
            chain_id = x[1]
            residue_num = x[2]
            #get the corresponding pdb file
            model = atomium.fetch(pdb_code).model
            chain = model.chain(chain_id)
            residue = chain.residue(f"{chain_id}.{residue_num}")

            #all the atoms in that residue of the feature to be compared
            ft_atoms = []
            for atom in residue.atoms():
               ft_atoms.append(atom.id)

            #now have list of atoms in that residue
            #get all the atoms in the residue of interest in this pdb/chain id

            for i in pdb_infos_res:
               #only get distances for the residues numbers from the same pdb file
               if i[0] == pdb_code and i[1] == chain_id:
                  the_residue = i[2]

                  #get a list of atoms in that residue
                  res_atoms = []
                  residue = chain.residue(f"{chain_id}.{the_residue}")
                  for atom in residue.atoms():
                     res_atoms.append(atom.id)

                  #get the distance between each atom combination
                  for a in ft_atoms:
                     for b in res_atoms:
                        d = model.atom(int(a)).distance_to(model.atom(int(b)))
                        #record the smallest distance
                        if d < min_dist:
                           min_dist = d
                           relevant_pdb = pdb_code
                           relevant_chain = chain_id
                           relevant_ft_residue = residue_num


      #record the smallest distance and the related info for that feature into the list 
      feature_infos_pdb.append([min_dist, relevant_pdb, relevant_chain, relevant_ft_residue])

   #loop through all the features in the pdb_infos_fts list
   #list with an item for ech feature
   return (feature_infos_pdb)


#*************************************************************************
def read_url_pdbsws (pdb, chain, pdb_residue):
   """ Get UniProt id and residue number equivalents from a specified PDB 
   residue.

   Input:   pdb            --- 4 letter PDB code
            chain          --- Chain ID
            pdb_residue    --- PDB Residue number
    Output:  file            --- Text file from url as a byt

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


