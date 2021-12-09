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
sprot_file = read_url_sprot(uniprot_ac)


#get list of residues in features
sp_ft_residues = get_ft_residues(sprot_file, res_of_interest)

#get pdb code and residue number for the protein and residue of interest

#list of PDB code, chain and residue number 
pdb_info_res = pdb_sws(uniprot_ac, res_of_interest)
#if doesnt work - bc res_of interest is an integer not a string from command line

#get pdb residue numbers for the feature residues



#get coordinates of atoms in residue of interest



#get coordinates of atoms for residues in features



#calculate distances between atoms in residue of interest and atoms in feature residues
