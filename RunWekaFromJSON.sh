#!/bin/sh

#*************************************************************************
#
#   Program:    RunWekaFromJSON
#   File:       RunWekaFromJSON.sh
#   
#   Version:    V1.0
#   Date:       15.03.22
#   Function:   Runs Weka Random Forest train and test from JSON files
#
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
#   Usage: RunWekaFromJSON.sh fold NTREE NFT DEPTH 
#           
#           fold    - number for x-fold cross validation
#           NTREE   - Machine Learning parameter tree
#           NFT     - Machine Learning parameter number of features
#           DEPTH   - Machine Learning paramtere depth
#
#   Required:
#       files:
#           inputs.dat/inputs_updated.dat
#       programs:
#           json2csv_uniprot_allPDB.pl
#           join_csv.sh
#           xvalidate.pl
#           CreateTrainTest.pl 
#           csv2arff
#           
#           weka.classifiers.trees.RandomForest
#           
#
#*************************************************************************
#
#   Revision History:
#   =================
#   V1.0   15.03.22  Original       By: BAM
#   V1.1   15.03.22  Fixed mistakes By: ACRM
#*************************************************************************


# cross validation number from command line
xval=$1

# Scripts
bindir=/home/bmikucka/bin

# Data
PD=/home/bmikucka/data/PD/
SNP=/home/bmikucka/data/SNP/

# Weka 
export WEKA=/home/amartin/weka-3-8-3
export CLASSPATH="$WEKA/weka.jar"
CLASSIFIER="weka.classifiers.trees.RandomForest"
JAVAFLAGS=-Xmx6g
NTREE=$2
NFEAT=$3
DEPTH=$4

# Storage of all Intermediate files
mkdir -p INTERMEDIATEFILES

echo -n "Converting json to csv for PDs..."
mkdir -p pdcsv
for file in $PD/*.json
do
   echo $file >> json2csv_pd.log
   $bindir/cleanjson.pl $file > cleanedp.json
   $bindir/json2csv_uniprot_allPDB.pl cleanedp.json PD >pdcsv/`basename $file .json`.csv 2>> json2csv_pd.log
   \rm cleanedp.json
done
echo "done"

echo -n "Converting json to csv for SNPs..."
mkdir -p snpcsv
for file in $SNP/*.json
do
   echo $file >> json2csv_snp.log
   $bindir/cleanjson.pl $file > cleaneds.json
   $bindir/json2csv_uniprot_allPDB.pl cleaneds.json SNP >snpcsv/`basename $file .json`.csv 2>> json2csv_snp.log
   \rm cleaneds.json
done
echo "done"

echo -n "Removing PD error files..."
for file in pdcsv/*
do
    if grep -q ERROR $file ; then
        rm -f $file
    fi
done
echo "done"

echo -n "Removing SNP error files..."
for file in snpcsv/*
do
    if grep -q ERROR $file ; then
        rm -f $file
    fi
done
echo "done"

echo -n "Concatenating csv files..."
$bindir/join_csv.sh pd 
$bindir/join_csv.sh snp
echo "done"

#echo -n "Editing SProtFT-dist columns >15A and -1 to 15A... "
#$bindir/fix_csv.py pd.csv
#$bindir/fix_csv.py snp.csv
#echo "done"

#find dataset size for PD
numPD=$(cat pd.csv | wc -l)
numPD=$(expr $numPD - 1)

# find dataset size for SNP
numSNP=$(cat snp.csv | wc -l)
numSNP=$(expr $numSNP - 1)

#number of balancing runs depending on how imbalanced the datasets are
if [ $numSNP -gt $numPD ];
then
   nBalance=$(expr $numSNP / $numPD)
else
   nBalance=$(expr $numPD / $numSNP)
fi

echo "Completed calculations of datapoints and number of balancing runs"
echo "$nBalance balancing runs will be performed"

#split the files into subdirectories
echo -n "Creating $xval cross-validation folds for PDs..."
mkdir -p PD_folds
cd PD_folds
$bindir/xvalidate.pl -n=$xval ../pd.csv 
cd ..
echo "done"

echo -n "Creating $xval cross-validation folds for SNPs..."
mkdir -p SNP_folds
cd SNP_folds
$bindir/xvalidate.pl -n=$xval ../snp.csv 
cd ..
echo "done"

echo -n "Combining PD and SNP folds..."
mkdir -p JOIN_folds
for file in PD_folds/*.csv
do
    csvfile=`basename $file`
    joinfile=JOIN_folds/$csvfile
    head -1 $file >$joinfile
    grep -v Binding PD_folds/$csvfile  >> $joinfile
    grep -v Binding SNP_folds/$csvfile >> $joinfile
done
echo "done"


echo -n "Creating $xval training and testing sets..."
for ((fold=0; fold<$xval; fold++))
do
    $bindir/CreateTrainTest.pl $fold $xval JOIN_folds
done
echo "done"

echo "Doing the training and testing runs!"

fold=0
for trainCSV in train_*.csv
do
    # Find the smaller dataset size to set the limit in training for the balancing runs
    numPD=$(grep PD $trainCSV | wc -l)
    numPD=$(expr $numPD - 1)
    numSNP=$(grep SNP $trainCSV | wc -l)
    numSNP=$(expr $numSNP - 1)

    if [ $numSNP -gt $numPD ];
    then
        limit=$numPD
    else
        limit=$numSNP
    fi

    # Find the name of the test file
    testCSV=`echo $trainCSV | sed s/train/test/`

    # Create the test ARFF file
    testARFF=`basename $testCSV .csv`.arff
    errors=`basename $testCSV .csv`.errors
    echo ""
    echo ""
    echo -n "FOLD $fold - Creating test ARFF file: $testARFF..."
    csv2arff -skip -ni inputs.dat dataset $testCSV > $testARFF 2>$errors
    echo "done"
    
    # Do the balancing training runs
    for ((balance=0; balance<$nBalance; balance++))
    do
        echo ""
        echo "RUNNING FOLD $fold BALANCING RUN $balance"

	# Convert training data to arff format
        trainARFF=`basename $trainCSV .csv`_$balance.arff
        errors=`basename $trainCSV .csv`_$balance.errors
        echo -n "Creating train ARFF file: $trainARFF..."
	csv2arff -skip -ni -limit=$limit inputs.dat dataset $trainCSV >$trainARFF 2>$errors
	# -skip      - skip records with missing values
	# -ni        - do not convert binary inputs to nominal Boolean
	# -limits=n  - balancing dataset 
	# inputs.dat - file containing list of input fields
	# dataset    - the name of the output field in the CSV file
	# train.csv  - the input csv file
	# train.arff - the output arff file
        echo "done"
        
	# This trains, saves the trained model and tests in one go
	# To train without the independent testing remove the '-T test.arff'
        trainModel=`basename $trainCSV .csv`_$balance.model
        echo -n "Training and testing..."
        outFile=test_${fold}_${balance}.out
	java $JAVAFLAGS $CLASSIFIER -I $NTREE -K $NFEAT -depth $DEPTH -t $trainARFF -d $trainModel -T $testARFF > $outFile
        echo "done (Results in $outFile)"

	# This tests on a pre-trained model
	#java $CLASSIFIER -l $trainModel -T $testARFF -p 0 > test2_${fold}_${balance}.out
    done
    fold=$(expr $fold + 1)
    mv *.arff *.errors *.model INTERMEDIATEFILES
done

