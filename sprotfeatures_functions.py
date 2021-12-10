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
   Output:  sp_ft_residues     --- List of SwissProt resdiue numbers in 
                                relevant features
           

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
                  #only the specific residues are involved
                  sp_ft_residues.append(start)
                  sp_ft_residues.append(stop)

                  
               elif len(info_list) == 4:
                  #for features with residue numbers separated by spaces
                  start = int(info_list[2])
                  stop = int(info_list[3])
                  #add the two residues to the sp_ft_residues list
                  sp_ft_residues.append(start)
                  sp_ft_residues.append(stop)


            else:
               #for features with range of residues
               if '..' in info_list[2]:
                  res_range_str = info_list[2].replace('..', ' ')
                  res_range = res_range_str.split() 
                  start = int(res_range[0])
                  stop = int(res_range[1])
                  #get numbers of all the residues
                  residues = range(start, stop)
                  #add everything from start to stop to the sp_ft_residues list
                  sp_ft_residues.extend(residues)
                  
               else:

                  #for features at one residue
                  if len(info_list) == 3:
                     start = int(info_list[2])
                     stop = int(info_list[2])
                     #add the residue to the sp_ft_residues list
                     sp_ft_residues.append(start)

                  #for features with residue numbers separated by spaces
                  elif len(info_list) == 4:
                     start = int(info_list[2])
                     stop = int(info_list[3])
                     #get numbers of all the residues
                     residues = range(start, stop)
                     #add everything from start to stop to the sp_ft_residues list
                     sp_ft_residues.extend(residues)

            
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
def read_url_pdbsws (uniprot_ac, uniprot_resid):
   """Reads text file from url PDBSWS to get information on PDB files 
   related to the UniProt accession number.

   Input:   uniprot_ac      --- UniProt accession number  
            uniprot_resid   --- Number for residue of interest in UniProt file
   Output:  file            --- Text file from url as a string
   
   09.12.21    Original    By: BAM

   """
  
   url = 'http://www.bioinf.org.uk/servers/pdbsws/query.cgi?plain=1&qtype=ac&id={}&res={}'.format(uniprot_ac, uniprot_resid)
   file = urlopen(url)
   return (file.read())


#*************************************************************************

def pdb_sws (uniprot_ac, uniprot_resid):
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

   pdb_file_byt = read_url_pdbsws(uniprot_ac, uniprot_resid)
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
def get_best_distance (pdb_infos_res, pdb_infos_fts):
   """ Processes lists of information from PDB about residue of interest 
   and features, gets shortest distance between a residue of interest atom
   and one of the feature atoms.

   Input:   pdb_infos_res --- list of lists with 3 items (PDB code, chain
                              identifier, and residue number) for the residue 
                              of interest
            pdb_infos_fts --- list of lists for each feature with lists ith 
                              3 items (PDB code, chain identifier, and residue 
                              number) for the features
   Output:  d             --- shortest distance between 2 atoms (Angstrom)
            ft_residue    --- the feature residue (PDB number) closest to 
                              residue of interest
            ft_atom       --- the feature residue atom closest to 
                              residue of interest
            res_atom      --- the atom number in the residue of interest 
                              closest to the closest feature
   
   10.12.21    Original    By: BAM

   """

   #unpack the list to get info on residue of interest
   for x in pdb_infos_res:
      pdb_code = x[0]
      chain_id = x[1]
      residue_num = x[2]
      #get the corresponding pdb file
      model = atomium.fetch(pdb_code).model
      #is this case sensitive?
      chain = model.chain(chain_id)
      residue = chain.residue(f"{chain_id}.{residue_num}")

      #make a list of feature residue numbers that will be compared in 
      #this PDB and chain combination
      ft_residues = []

      #make a lits of atoms corresponding to that residue 
      res_atoms = []
      for atom in residue.atoms():
         res_atoms.append(atom.id)

      for y in pdb_infos_fts:
         #for each residue from a feature
         
         for z in y:
            #for each list with info about that residue (list in a list)
            if z[0] == pdb_code and z[1] == chain_id:
               #compare only between the same PDBs and chains
               #add the residue number into the list
               ft_residues.append(z[2])

      #now have feature residue number list and the residue of interest number

      #dictionary with residue numbers as keys and atom numbers in lists as values
      ft_atoms = {}
      atoms = []
      #get atom number for each residue
      for i in ft_residues:
         residue = chain.residue(f"{chain_id}.{i}")
         for atom in residue.atoms():
            atoms.append(atom.id)

         #add this to the dictionary
         ft_atoms[i] = atoms

      d = 50 #start distance
      #could have this as a variable

      #for this residue number get distance between all its atoms and 
      #atoms of feature residues
      for i in res_atoms:
         #i is an atom of residue of interest
         i_int = int(i)
         for key in ft_atoms:
            #key is a residue in features
            for j in ft_atoms[key]:
               #ft_atoms[key] is the list of atoms - j is one atom number
               #get the distance between atom of interest and atom in feature residue
               j_int = int(j)
               dd = pdb_code.model.atom(i_int).distance_to(pdb_code.model.atom(j_int))
               if dd < d:
                  #if this distance is smaller than any of the previous ones recorded:
                  d = dd 
                  ft_atom = j
                  ft_residue = key
                  res_atom = i
                  relevant_pdb = pdb_code
                  relevant_chain = chain_num


   return (d, ft_residue, ft_atom, res_atom, relevant_chain, relevant_pdb)


#*************************************************************************
