

"""
Program: sprotfeatures
File:    sprotfeatures.py

Version:    V1.0
Date:       09.12.2021
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
from urllib.request import urlopen
import atomium
from sprotfeatures_functions import (read_file, read_url_sprot, 
   get_ft_residues, check_feature, read_url_pdbsws, pdb_sws, get_best_distance)

#*************************************************************************

#take uniprot accession number and residue number from command line

print('running the program...')

uniprot_ac = sys.argv[1]
res_of_interest = int(sys.argv[2])


#get SwissProt file based on uniprot_ac
sprot_file_byt = read_url_sprot(uniprot_ac)
#this is now a byte object - convert into string
encoding = 'utf-8'
sprot_str = sprot_file_byt.decode(encoding)

#run check_feature first?

#get list of residues in features
#res_of_interest an integer here
sp_ft_residues = get_ft_residues(sprot_str, res_of_interest)


#get pdb code and residue number for the protein and residue of interest
#list of PDB codes, chains and residue numbers 
pdb_infos_res = pdb_sws(uniprot_ac, res_of_interest)
#works for P03952(uniprot - works for both) but not for Q6GZV6(swissprot-only for sprot)

#remove repeats 
temp_list = []
for i in pdb_infos_res:
    if i not in temp_list:
        temp_list.append(i)
pdb_infos_res = temp_list


#get pdb residue numbers for the feature residues
pdb_infos_fts = []
for residue in sp_ft_residues:
   #list of PDB infos for one residue
   pdb_infos_ft = pdb_sws(uniprot_ac, residue)
   #combine lists for a PDB info for all feature residues list
   pdb_infos_fts.append(pdb_infos_ft)
   #will have empty lists if the SwissProt residue doesn't have a PDB residue equivalent

#remove repeats
temp_list = []
for i in pdb_infos_fts:

   if i not in temp_list:
      temp_list.append(i)
   #remove empty lists (when SwissProt residues are not in the PDB file)
   if i == [['', '', '']]:
      temp_list.remove(i)

pdb_infos_fts = temp_list


#get the shortest distance between an atom in residue of interest and
#an atom in a residue in a feature 

(d, ft_residue, ft_atom, res_atom, relevant_chain, relevant_pdb) = get_best_distance (pdb_infos_res, pdb_infos_fts)

print (d, ft_residue, ft_atom, res_atom, relevant_chain, relevant_pdb)

#d: shortest distance between atoms
#ft_residue: PDB residue number of the feature residue closest to the residue of interest
#ft_atom: atom number 
#res_atom: atom in the residue of interest that is closest to the feature


#from ft_residue (PDB numbers) get what the feature is
#from PBD to SwissProt

#From SwissProt to feature name

print ('done running')



