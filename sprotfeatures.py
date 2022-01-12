#!/usr/bin/python3

"""
Program: sprotfeatures
File:    sprotfeatures.py

Version:    V1.0
Date:       09.12.2021
Function:   Calculates the distance of mutant residue to the closest 
            SwissProt feature. 


Author: Barbara A. Mikucka

--------------------------------------------------------------------------
Description:


--------------------------------------------------------------------------
Usage:
sprotfeatures.py [chain]resnum[insert] newaa pdbfile


--------------------------------------------------------------------------
Revision history:
=================
V1.0  16.11.21    Original    By: BAM
"""

#*************************************************************************
# Import libraries

import sys
import functools 
import timeit
import re
from urllib.request import urlopen
import atomium
import json
from sprotfeatures_functions import (read_file, process_resnum, get_pdb_code, read_url_sprot, 
   get_ft_residues, check_feature, read_url_swspdb, sws_pdb, feature_distance, 
   read_url_pdbsws, pdb_sws, residue_to_feature)

#*************************************************************************

#take uniprot accession number and residue number from command line

res_id = [1]
newaa = [2]
pdbfile = [3]

#uniprot_ac = sys.argv[1]
#mut_res_number = sys.argv[2]
#mutant_residue_id = sys.argv[3]

#uniprot_ac, mut_res_number, mutant_residue_id

@functools.lru_cache (maxsize = None)
def get_results (res_id, newaa, pdbfile):

   #get PDB code
   pdb_code = get_pdb_code (pdbfile)

   #get the residue number, chain and insert from resnum
   (chain_id, resnum_pdb, insert) = process_resnum (res_id)

   #get SwissProt accession number 
   (sprot_ac, resnum_sprot) = pdb_sws (pdb_code, chain_id, resnum_pdb)

   #get SwissProt file based on uniprot_ac
   sprot_file_byt = read_url_sprot(sprot_ac)
   #this is now a byte object - convert into string
   encoding = 'utf-8'
   sprot_str = sprot_file_byt.decode(encoding)

   #get list of residues in features
   #mut_res_number an integer here
   sp_ft_residues = get_ft_residues(sprot_str, int(resnum_sprot))
   #this is a list of lists - one list per feature


   #get pdb code and residue number for the protein and residue of interest
   #list of PDB codes, chains and residue numbers 
   pdb_infos_res = sws_pdb(sprot_ac, resnum_sprot)

   #remove repeats 
   temp_list = []
   for i in pdb_infos_res:
       if i not in temp_list:
           temp_list.append(i)
   pdb_infos_res = temp_list


   #get pdb residue numbers for the feature residues
   pdb_infos_fts = []
   for feature in sp_ft_residues:
      #for each list of a few residues in the feature
      per_feature = []
      
      for residue in feature: 
         #for each of the residues in that feature
         per_residue = []
         #list of PDB infos for one residue
         per_residue = sws_pdb (uniprot_ac, residue)
         #from the function, a list of lists is added that corresponds to that residue number

         if per_residue != []:
            #combine lists for a PDB info for all residue in one feature
            per_feature.append(per_residue)

      #combine all non-empty into a list for all the features
      if per_feature != []:
         pdb_infos_fts.append(per_feature)



   #get list with the distances for each feature with the relevant information
   feature_infos_pdb = feature_distance (pdb_infos_res, pdb_infos_fts)

   #dictionary with all the feature names as keys and corresponding 
   #distances from the residue of interest 
   final_infos = {}

   for feature in feature_infos_pdb:
      #information on where the shortest distance is between residue and feature
      min_dist = round(float(feature [0]), 3)
      relevant_pdb = feature [1]
      relevant_chain = feature [2]
      relevant_ft_residue = feature [3]

      #get the SwissProt information
      (sprot_ac, sprot_residue) = pdb_sws (relevant_pdb, relevant_chain, relevant_ft_residue)

      #get feature name for this feature
      feature_id = residue_to_feature (sprot_residue, sprot_ac)

      #add to the list of the info including the feature name if this is 
      #shorter than the one that's already there 
      if feature_id in final_infos:
         if min_dist < final_infos[feature_id]:
            final_infos[feature_id] = min_dist
      else:
         final_infos[feature_id] = min_dist

   #convert into json format
   output = json.dumps(final_infos)

   return output

output = get_results (uniprot_ac, mut_res_number, mutant_residue_id)
#print(timeit.timeit ('get_results (uniprot_ac, mut_res_number, mutant_residue_id)', globals=globals(), number=1))

print(output)




