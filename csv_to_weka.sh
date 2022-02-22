
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

lines_pd=$(expr $N / $x)
lines_snp=$(expr $n / $x)

# change these to:
# divide into x folds (x from command line) - make sure no overlap in proteins between folds 
# combine snp and pd 

# split into x parts 
# PD
gsplit -a 4 -d -l $(lines_pd) pd.csv pd_part.csv 
#get rid of first line in first file
cat pd_part.csv0000 > tmp.csv
grep -v Binding tmp.csv >> pd_part.csv0000 


# check if there is a file with leftovers - add to last file with if there
file_number=$(echo ${x} | awk '{ printf "%04i\n", $0 }')
file_number_2=$(echo $(expr(${x} - 1)) | awk '{ printf "%04i\n", $0 }')

if test -f "pd_part.csv${file_number}"; then
	cat "pd_part.csv${file_number}" >> "pd_part.csv${file_number_2}"
fi

# check the same protein is not in different groups
for i in {1..$x}
	# get names of two consecutive files
	file_num_2=$(echo ${i} | awk '{ printf "%04i\n", $0 }')
	file_num_1=$(echo $(expr(${i} - 1)) | awk '{ printf "%04i\n", $0 }')
	file_1="pd_part.csv${file_num_1}"
	file_2="pd_part.csv${file_num_2}"

	for (( a=0; a=1; ))
	do  
		pdb_last=$(awk 'END{print substr($0,0,4)}' $file_num_1)
		pdb_first=$(cut -c 1-4 $file_num_2)
		# if the PDB code in the last line of file 1 and the first line 
		# of file 2 are the same then move the last line to file 2.
		if $pdb_last==$pdb_first; 
			then tail -n 1 "${file_1}" >> "${file_2}"
			else a=1
		fi
	done
done

# SNP
gsplit -a 4 -d -l $(lines_snp) snp.csv snp_part.csv 
#get rid of first line in first file
cat snp_part.csv0000 > tmp.csv
grep -v Binding tmp.csv >> snp_part.csv0000 

# check if there is a file with leftovers
if test -f "snp_part.csv${file_number}"; then
	cat "snp_part.csv${file_number}" >> "snp_part.csv${file_number_2}"
fi

# check the same protein is not in different groups
for i in {1..$x}
	# get names of two consecutive files
	file_num_2=$(echo ${i} | awk '{ printf "%04i\n", $0 }')
	file_num_1=$(echo $(expr(${i} - 1)) | awk '{ printf "%04i\n", $0 }')
	file_1="snp_part.csv${file_num_1}"
	file_2="snp_part.csv${file_num_2}"

	for (( a=0; a=1; ))
	do  
		pdb_last=$(awk 'END{print substr($0,0,4)}' $file_num_1)
		pdb_first=$(cut -c 1-4 $file_num_2)
		# if the PDB code in the last line of file 1 and the first 
		# line of file 2 are the same then move the last line to file 2.
		if $pdb_last==$pdb_first; 
			then tail -n 1 "${file_1}" >> "${file_2}"
			else a=1
		fi
	done
done


# join pd and snp with same suffix
for i in pd_part.csv*; do
   p="${i#pd_}"
   [[ -f "snp_$p" ]] && cat "$i" "snp_$p" > "$p"
done
# now have files with snp and pd data - part.csv0000 (x number)


for i in {1..$x}
do
	head -1 pd.csv > train.csv 
	head -1 pd.csv > test.csv

	#join all files but leave one testing set
	test_num=$(echo ${i} | awk '{ printf "%04i\n", $0 }')
	test_file="part.csv${test_num}"

	# make test and train files
	cat test_file >> tmp_test.csv
	cat part.csv* !($test_file) > tmp_train.csv

	#clean up the files
	grep -v Binding tmp_train.csv >> train.csv
	grep -v Binding tmp_test.csv >> test.csv

	#repeat this step multiple times - balancing

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
	java $CLASSIFIER -I $NTREE -t train.arff -d train.model -T test.arff > test_${i}.out

	# This tests on a pre-trained model
	java $CLASSIFIER -l train.model -T test.arff > test2_${i}.out


	# Cleanup
	rm train.csv test.csv train.arff test.arff

done

