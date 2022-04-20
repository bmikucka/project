 #!/bin/sh

#*************************************************************************
#
#   Program:    join_csv
#   File:       join_csv.sh
#   
#   Version:    V1.0
#   Date:       23.02.22
#   Function:   Join CSV files
#
#   Author: Barbara A. Mikucka
#               
#*************************************************************************
#
#   Description:
#   ============
#
#*************************************************************************
#
#   Usage:
#   ======
#
#*************************************************************************
#
#   Revision History:
#   =================
#   V1.0   23.02.22  Original   By: BAM
#
#*************************************************************************

#join all csvs together - only keep one header
#take dataset type SNP/PD from command line

dataset=$1

# Join the csv files
   cat ${dataset}csv/* > tmp.csv
   head -1 tmp.csv > ${dataset}.csv
   grep -v Binding tmp.csv >> ${dataset}.csv
   rm tmp.csv
