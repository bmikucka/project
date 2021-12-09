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
from sprotfeatures_functions import (read_file, check_feature, pdb_sws)

#*************************************************************************

#take file and residue number from command line

#sprot_file = sys.argv[1]
#res_of_interest = int(sys.argv[2])

#print (check_feature (sprot_file, res_of_interest))

pdb_sws(P03952, 400)


