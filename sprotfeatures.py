#!/usr/bin/python3

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
from urllib.request import urlopen
from sprotfeatures_functions import (read_file, read_url_sprot, 
   get_ft_residues, check_feature, read_url_pdbsws, pdb_sws)

#*************************************************************************

#take uniprot accession number and residue number from command line

uniprot_ac = sys.argv[1]
res_of_interest = int(sys.argv[2])


#get SwissProt file based on uniprot_ac
sprot_file_byt = read_url_sprot(uniprot_ac)
#this is now a byte object - convert into string
encoding = 'utf-8'
sprot_str = sprot_file_byt.decode(encoding)


#get list of residues in features
#res_of_interest an integer here
sp_ft_residues = get_ft_residues(sprot_str, res_of_interest)


#get pdb code and residue number for the protein and residue of interest
#list of PDB codes, chains and residue numbers 
pdb_infos_res = pdb_sws(uniprot_ac, res_of_interest)
#works for P03952(uniprot - works for both) but not for Q6GZV6(swissprot-only for sprot)

#get pdb residue numbers for the feature residues
pdb_infos_fts = []
for residue in sp_ft_residues:
   #list of PDB infos for one residue
   pdb_infos_ft = pdb_sws(uniprot_ac, residue)
   #combine lists for a PDB info for all feature residues list
   pdb_infos_fts.append(pdb_infos_ft)

#will have empty lists if the SwissProt residue doesn't have a PDB residue equivalent


#get coordinates of atoms in residue of interest



#get coordinates of atoms for residues in features



#calculate distances between atoms in residue of interest and atoms in feature residues
