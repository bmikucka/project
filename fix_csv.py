#!/usr/bin/python3

"""
Program: fix_csv
File:    fix_csv.py

Version:    V1.0
Date:       15.03.2022
Function:   Changes all values in SprotFTdist columns in csv file to 15 if 
            value is >15.


Author: Barbara A. Mikucka

--------------------------------------------------------------------------
Usage:
=================


--------------------------------------------------------------------------
Revision history:
=================
V1.0  15.03.22    Original    By: BAM
"""

#*************************************************************************
# Import libraries

import sys
import pandas as pd

#*************************************************************************

filename = sys.argv[1]

#read csv file
df = pd.read_csv(filename)

#get number of rows
row_count = len(df)

# for column 1 get the value 
for x in range (row_count):
   value = df._get_value(x, 'SprotFTdist-ACT_SITE')
   if value != '?':
      distance = float(value)

      if distance > 15 or distance == -1:
         #change the value to 15
         df.loc[x, 'SprotFTdist-ACT_SITE'] = '15'
         # writing into the file
         df.to_csv(filename, index=False)

      distance = float(df._get_value(x, 'SprotFTdist-BINDING'))
      if distance > 15 or distance == -1:
         #change the value to 15
         df.loc[x, 'SprotFTdist-BINDING'] = '15'
         # writing into the file
         df.to_csv(filename, index=False)

      distance = float(df._get_value(x, 'SprotFTdist-CA_BIND'))
      if distance > 15 or distance == -1:
         #change the value to 15
         df.loc[x, 'SprotFTdist-CA_BIND'] = '15'
         # writing into the file
         df.to_csv(filename, index=False)

      distance = float(df._get_value(x, 'SprotFTdist-DNA_BIND'))
      if distance > 15 or distance == -1:
         #change the value to 15
         df.loc[x, 'SprotFTdist-DNA_BIND'] = '15'
         # writing into the file
         df.to_csv(filename, index=False)
     
      distance = float(df._get_value(x, 'SprotFTdist-NP_BIND'))
      if distance > 15 or distance == -1:
         #change the value to 15
         df.loc[x, 'SprotFTdist-NP_BIND'] = '15'
         # writing into the file
         df.to_csv(filename, index=False)

      distance = float(df._get_value(x, 'SprotFTdist-METAL'))
      if distance > 15 or distance == -1:
         #change the value to 15
         df.loc[x, 'SprotFTdist-METAL'] = '15'
         # writing into the file
         df.to_csv(filename, index=False)

      distance = float(df._get_value(x, 'SprotFTdist-MOD_RES'))
      if distance > 15 or distance == -1:
         #change the value to 15
         df.loc[x, 'SprotFTdist-MOD_RES'] = '15'
         # writing into the file
         df.to_csv(filename, index=False)

      distance = float(df._get_value(x, 'SprotFTdist-CARBOHYD'))
      if distance > 15 or distance == -1:
         #change the value to 15
         df.loc[x, 'SprotFTdist-CARBOHYD'] = '15'
         # writing into the file
         df.to_csv(filename, index=False)

      distance = float(df._get_value(x, 'SprotFTdist-MOTIF'))
      if distance > 15 or distance == -1:
         #change the value to 15
         df.loc[x, 'SprotFTdist-MOTIF'] = '15'
         # writing into the file
         df.to_csv(filename, index=False)

      distance = float(df._get_value(x, 'SprotFTdist-LIPID'))
      if distance > 15 or distance == -1:
         #change the value to 15
         df.loc[x, 'SprotFTdist-LIPID'] = '15'
         # writing into the file
         df.to_csv(filename, index=False)



