
#!/bin/sh

#needed to run:
#SNP and PD directories with relevant json files 
#join_csv.sh, inputs_updated.dat

# Input: a number x for x-fold validation

# Weka 
export WEKA=/home/amartin/weka-3-8-3
export CLASSPATH="$WEKA/weka.jar"

#convert json to csv
mkdir -p snpcsv
   for file in SNP/*.json
   do
   	   echo $file >> json2csv.log
   	   cleanjson.pl $file > cleaned.json
       ~/bin/json2csv_uniprot_allPDB_pl.txt $file SNP >snpcsv/`basename $file .json`.csv 2>> json2csv.log
       rm cleaned.json
   done

mkdir -p pdcsv
   for file in PD/*.json
   do
   	   echo $file >> json2csv.log
   	   cleanjson.pl $file > cleaned.json
       ~/bin/json2csv_uniprot_allPDB_pl.txt $file PD >pdcsv/`basename $file .json`.csv 2>> json2csv.log
       rm cleaned.json
   done

echo Done: JSON to CSV conversion.

#remove error files
for file in snpcsv/*
   do
       if grep -q ERROR $file ; then
           rm -f $file
       fi
   done

for file in pdcsv/*
   do
       if grep -q ERROR $file ; then
           rm -f $file
       fi
   done

echo Done: removed csv files with errors from pdcsv and snpcsv.


# join csv files 
~/bin/join_csv.sh snp
~/bin/join_csv.sh pd 

echo Done: All SNP and PD csv files joined into snp.csv and pd.csv.

# cross validation number from command line
x=$1

# find n - dataset size for SNP
SNP_lines=$(cat snp.csv | wc -l)
snp_datasize=$(expr $SNP_lines - 1)
#find N - dataset size for PD
PD_lines=$(cat pd.csv | wc -l)
pd_datasize=$(expr $PD_lines - 1)

#number of balancing runs depending on how imbalanced the datasets are
if [ $pd_datasize -gt $snp_datasize ];
	then
		balancing_runs=$(expr $pd_datasize / $snp_datasize)
		small_dataset=$pd_datasize
	else
		balancing_runs=$(expr $snp_datasize / $pd_datasize)
		small_dataset=$snp_datasize
fi

echo Done: Calculations of datapoints and number of balancing runs.

#split the files into subdirectories
mkdir -p PD_folds
cd PD_folds
~/bin/xvalidate.pl -n=${x} ../pd.csv 
cd ..


mkdir -p SNP_folds
cd SNP_folds
~/bin/xvalidate.pl -n=${x} ../snp.csv 
cd ..

echo Done: SNP and PD split into $x folds.

# join pd and snp with same suffix
mkdir -p JOIN_folds
ls PD_folds | while read file; do
  cat PD_folds/"$file" SNP_folds/"$file" >> JOIN_folds/tmp.csv
  head -1 JOIN_folds/tmp.csv > JOIN_folds/"$file"
  grep -v Binding JOIN_folds/tmp.csv >> JOIN_folds/"$file"
  rm JOIN_folds/tmp.csv
done

echo Done: SNP and PD folds joined into $x files- in JOIN_folds directory.

#combine folds with a different fold as test set each time
for ((i=0; i<$x; i++)); do

	head -1 pd.csv > train.csv 
	head -1 pd.csv > test.csv

	#join all files but leave one testing set
	test_num=$(expr $i - 1)
	#test_file= /*${test_num}.csv

	# make test and train files
	cat JOIN_folds/*${test_num}.csv > ./tmp_test.csv
	#cat JOIN_folds/fold_* !(JOIN_folds/*${test_num}.csv) > ./tmp_train.csv
	cat ./JOIN_folds/!(*${test_num}.csv) > ./tmp_train.csv


	#clean up the files
	grep -v Binding tmp_train.csv >> train.csv
	grep -v Binding tmp_test.csv >> test.csv

	#clean up
	rm tmp_test.csv tmp_train.csv

	echo Done: Test and training sets set up for Weka ($i out of $x).

	# repeat the next step multiple times - balancing done N/n times
	# Note that CSV2arff will do the balancing for you by using -limit

	csv2arff -skip -ni -limit=${small_dataset} inputs_updated.dat dataset test.csv >test.arff

	for ((j=0; j<$balancing_runs; j++)); do

		# Convert training data to arff format
		csv2arff -skip -ni -limit=${small_dataset} inputs_updated.dat dataset train.csv >train.arff
		# -skip      - skip records with missing values
		# -ni        - do not convert binary inputs to nominal Boolean
		# -limits=n  - balancing dataset 
		# inputs.dat - file containing list of input fields
		# dataset    - the name of the output field in the CSV file
		# train.csv  - the input csv file
		# train.arff - the output arff file
		# Convert test data to arff format

		# Set parameters for the machine learning
		CLASSIFIER="weka.classifiers.trees.RandomForest"
		NTREE=100

		# This trains, saves the trained model and tests in one go
		# To train without the independent testing remove the '-T test.arff'
		java $CLASSIFIER -I $NTREE -t train.arff -d train.model -T test.arff > test_${i}_${j}.out

		# This tests on a pre-trained model
		java $CLASSIFIER -l train.model -T test.arff > test2_${i}_${j}.out


		# Cleanup
		rm train.csv test.csv train.arff test.arff

		echo Done: Weka for training set $i and $j (cross-validation).
	done
done

