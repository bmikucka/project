#join all csvs together - only keep one header
#take SNP/PD from command line

dataset=$1

# Join the csv files
   cat ${dataset}csv/* > tmp.csv
   head -1 tmp.csv > ${dataset}.csv
   grep -v Binding tmp.csv >> ${dataset}.csv
