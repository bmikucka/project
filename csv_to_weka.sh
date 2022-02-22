
# Input: x-fold validation


# Note that CSV2arff will do the balancing for you by using -limit

# Weka 
export WEKA=/home/amartin/weka-3-8-3
export CLASSPATH="$WEKA/weka.jar"

# Split up the demo csv file into train and test sets

#cross validation number from command line
x=$1
#find n
number_1=$(wc -l snp.csv)
n=$(expr $number_1 - 1)
number_2=$(wc -l pd.csv)
N=$(expr $number_2 - 1)

# change these to:
# divide into x folds (x from command line) - make sure no overlap in proteins between folds 
# combine snp and pd 

# split into x parts 
# PD
gsplit -a 4 -d -l $(lines_pd) pd.csv pd_part.csv 
#check if there is a file with leftovers


# SNP
gsplit -a 4 -d -l $(lines_snp) snp.csv snp_part.csv 
#check if there is a file with leftovers


#join pd and snp with same suffix

for i in pd_part.csv*; do
   p="${i#pd_}"
   [[ -f "snp_$p" ]] && cat "$i" "snp_$p" > "$p"
done

#have files with snp and pd data - part.csv0000 (x number)

for i in {1..$x}
do
 
done

a=1
lines_pd=$(expr $N / $x)
lines_snp=$(expr $n / $x)




for i in {1..$x}
do
   head -1 pd.csv > pd${x}.csv
   
   >> pd${x}.csv


   head -1 snp.csv > snp${x}.csv

   a=$(expr $a + 1)
done



#check if (x+1)th file - if unequal split into x files
# make sure all but the first line are not headers


tail -11205 demo.csv >> test.csv

# I have created 'inputs.dat' to list the fields of interest as inputs
# for the machine learning

# Convert training data to arff format
csv2arff -skip -ni -limit=${n} inputs_updated.dat dataset train.csv >train.arff
# -skip      - skip records with missing values
# -ni        - do not convert binary inputs to nominal Boolean
# -limits=n  - balancing dataset 
# inputs.dat - file containing list of input fields
# dataset    - the name of the output field in the CSV file
# train.csv  - the input csv file
# train.arff - the output arff file
# Convert test data to arff format

csv2arff -skip -ni -limit=${n} inputs_updated.dat dataset test.csv >test.arff

# Set parameters for the machine learning
CLASSIFIER="weka.classifiers.trees.RandomForest"
NTREE=100

# This trains, saves the trained model and tests in one go
# To train without the independent testing remove the '-T test.arff'
java $CLASSIFIER -I $NTREE -t train.arff -d train.model -T test.arff > test.out

# This tests on a pre-trained model
java $CLASSIFIER -l train.model -T test.arff > test2.out


# Cleanup
rm train.csv test.csv train.arff test.arff