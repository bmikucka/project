#!/usr/bin/python3

"""
Program: sprotFTdist
File:    sprotFTdist.py

Version:    V1.0
Date:       09.12.2021
Function:   Calculates the distance of mutant residue to the closest 
            SwissProt feature. 


Author: Barbara A. Mikucka

--------------------------------------------------------------------------
Usage: sprotfeatures.py [chain]resnum[insert] newaa pdbfile [-vv] 
         [-nocache] [-force] [-info]

       (newaa maybe 3-letter or 1-letter code)
       -vv      Verbose
       -force   Force calculation even if results are cached
       -nocache Do not cache results
       -info    Prints a 1-line summary of what the plugin is doing


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
from sprotFTdist_lib import (os, read_file, process_resnum, 
   get_pdb_code, read_url_sprot, get_sprot_str, get_ft_residues, 
   pdb_ft_list, check_feature, read_url_swspdb, sws_pdb, feature_distance, 
   read_url_pdbsws, pdb_sws, residue_to_feature, get_bool_results,
   check_cache, write_cache)

#*************************************************************************

#take uniprot accession number and residue number from command line

opts = [opt for opt in sys.argv[1:] if opt.startswith("-")]
args = [arg for arg in sys.argv[1:] if not arg.startswith("-")]

res_id = args[0]
newaa = args[1]
pdbfile = args[2]

if '-info' in opts:
   print ("Calculating distance the mutant residue to feature residues")
   sys.exit()

if '-h' in opts:
   print ('''
   sprotFTdist.py, V1.0, 25.01.2022, Barbara A. Mikucka
   Function:   SAAPdap plugin: Calculates the distance of mutant residue to 
               the closest SwissProt feature. 

   --------------------------------------------------------------------------
   Usage: sprotfeatures.py [chain]resnum[insert] newaa pdbfile [-vv] 
   [-nocache] [-force] [-info]

      (newaa maybe 3-letter or 1-letter code)
      -vv      Verbose
      -force   Force calculation even if results are cached
      -nocache Do not cache results
      -info    Prints a 1-line summary of what the plugin is doing
            ''')
   sys.exit()


cached_file = check_cache (res_id, newaa, pdbfile, opts)
if cached_file != '':
   #print the results in json format
   print (cached_file)
elif cached_file == '':

   #get PDB code
   pdb_code = get_pdb_code (pdbfile)

   #get the residue number, chain and insert from resnum
   (chain_id, resnum_pdb, insert) = process_resnum (res_id)

   #get SwissProt accession number 
   (sprot_ac, resnum_sprot) = pdb_sws (pdb_code, chain_id, resnum_pdb)

   #get SwissProt file based on uniprot_ac
   sprot_str = get_sprot_str (sprot_ac)

   #get list of residues in features
   sp_ft_residues = get_ft_residues(sprot_str, int(resnum_sprot))
   #this is a list of lists - one list per feature

   #translate SwissProt numbering to PDB numbering
   pdb_ft_residues = pdb_ft_list (sp_ft_residues, sprot_ac, pdb_code, chain_id)

   #get list with the distances for each feature with the relevant information
   feature_distances = feature_distance (pdb_code, chain_id, resnum_pdb, pdb_ft_residues)


   #dictionary with all the feature names as keys and corresponding 
   #distances from the residue of interest 
   output = {
      "SprotFTdist-BOOL": 'tbd',
      'ACT_SITE':    -1,
      'BINDING':     -1, 
      'CA_BIND':     -1, 
      'DNA_BIND':    -1, 
      'NP_BIND':     -1, 
      'METAL':       -1, 
      'MOD_RES':     -1, 
      'CARBOHYD':    -1, 
      'MOTIF':       -1, 
      'LIPID':       -1
      }


   for feature in feature_distances:
      #information on where the shortest distance is between residue and feature
      min_dist = round(float(feature [0]), 3)
      closest_res_pdb = feature [1]

      #get the SwissProt information
      (sprot_ac, closest_res_sp) = pdb_sws (pdb_code, chain_id, closest_res_pdb)

      #get feature name for this feature
      feature_id = residue_to_feature (closest_res_sp, sprot_ac)

      #change the value for the feature if there isn't one there (=-1)
      if output[feature_id] != -1:
         #only replace if the distance is smaller than what is entered
         if min_dist < output[feature_id]:
            output[feature_id] = min_dist
      else:
         #for features that do not have a distance value yet
         output[feature_id] = min_dist

   distances_list = list(output.values())
   #exclude the BOOL result from distance evaluations
   distances_list = distances_list [1:]
   bool_result = get_bool_results(distances_list)

   #edit the BOOL to include the result
   output["SprotFTdist-BOOL"] = bool_result

   #change key names to fit output
   output['SprotFTdist-ACT_SITE'] = output.pop('ACT_SITE')
   output['SprotFTdist-BINDING'] = output.pop('BINDING')
   output['SprotFTdist-CA_BIND'] = output.pop('CA_BIND')
   output['SprotFTdist-DNA_BIND'] = output.pop('DNA_BIND')
   output['SprotFTdist-NP_BIND'] = output.pop('NP_BIND')
   output['SprotFTdist-METAL'] = output.pop('METAL')
   output['SprotFTdist-MOD_RES'] = output.pop('MOD_RES')
   output['SprotFTdist-CARBOHYD'] = output.pop('CARBOHYD')
   output['SprotFTdist-MOTIF'] = output.pop('MOTIF')
   output['SprotFTdist-LIPID'] = output.pop('LIPID')

   #convert into json format
   output = json.dumps(output)

   #cache results
   if '-nocache' not in opts:
      write_cache (res_id, newaa, pdbfile, output)

   print(output)
