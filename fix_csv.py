#!/usr/bin/python3

import sys
import pandas as pd

filename = sys.argv[1]

df = pd.read_csv(filename)

#get number of rows
row_count = len(df)

# for column 1 get the value 
for x in range (row_count):
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



