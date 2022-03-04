
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
       ~/bin/json2csv_uniprot_allPDB_pl.txt $file >snpcsv/`basename $file .json`.csv
   done

mkdir -p pdcsv
   for file in PD/*.json
   do
       ~/bin/json2csv_uniprot_allPDB_pl.txt $file >pdcsv/`basename $file .json`.csv
   done

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

# files ordered by SwissProt accession number 


# join csv files 
~/bin/join_csv.sh snp
~/bin/join_csv.sh pd 


# cross validation number from command line
x=$1

# find n - dataset size for SNP
number_1=$(wc -l snp.csv)
n=$(expr $number_1 - 1)
#find N - dataset size for PD
number_2=$(wc -l pd.csv)
N=$(expr $number_2 - 1)

#number of balancing runs depending on how imbalanced the datasets are
m=$(expr $N / $n)


#split the files into subdirectories
mkdir -p PD_folds
cd PD_folds
~/bin/xvalidate.pl -x${x} ../pd.csv 
cd ..


mkdir -p SNP_folds
cd SNP_folds
~/bin/xvalidate.pl -x${x} ../snp.csv 
cd ..

# join pd and snp with same suffix
mkdir -p JOIN_folds
ls PD_folds | while read file; do
  cat PD_folds/"$file" SNP_folds/"$file" >> JOIN_folds/tmp.csv
  head -1 JOIN_folds/tmp.csv > JOIN_folds/"$file"
  grep -v Binding JOIN_folds/tmp.csv >> JOIN_folds/"$file"
  rm JOIN_folds/tmp.csv
done


#combine folds with a different fold as test set each time
for i in {1..$x}
do
	head -1 pd.csv > train.csv 
	head -1 pd.csv > test.csv

	#join all files but leave one testing set
	test_num=$(expr ${i} - 1)
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

	# repeat the next step multiple times - balancing done N/n times
	# Note that CSV2arff will do the balancing for you by using -limit

	for j in {1..$m}
	do
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
		java $CLASSIFIER -I $NTREE -t train.arff -d train.model -T test.arff > test_${i}_${j}.out

		# This tests on a pre-trained model
		java $CLASSIFIER -l train.model -T test.arff > test2_${i}_${j}.out


		# Cleanup
		rm train.csv test.csv train.arff test.arff
	done
done

