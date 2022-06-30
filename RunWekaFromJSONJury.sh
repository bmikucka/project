#!/bin/sh

# Scripts
bindir=/serv/data/people/bmikucka/bin


# Weka 
export WEKA=/home/amartin/weka-3-8-3
export CLASSPATH="$WEKA/weka.jar"
CLASSIFIER="weka.classifiers.trees.RandomForest"
JAVAFLAGS=-Xmx6g



echo "Doing the testing runs!"


for ((fold=0; fold<10; fold++))
do
    echo "Handling fold no. ${fold}"

    #name of test file
    testCSV="test_0${fold}.csv"

    echo "Making individual csv files for test_0${fold}.csv datapoints"
    #divide CSV into individual CSVs
    mkdir test_files_${fold}
    $bindir/break_csv.sh $testCSV $fold

    echo "Creating ARFF files for all test files."

    #for each datapoint in the testing set - for each file test_0${fold}_${element}.csv
    element=0
    for file in ./test_0${fold}_*.csv
    do
        
        echo "Testing element ${element}"

        singletestCSV="test_0${fold}_${element}.csv"
        testARFF="test_0${fold}_${element}.arff"
        errors="test_0${fold}_${element}.errors"
        
        csv2arff -skip -ni inputs.dat dataset $singletestCSV >$testARFF 2>$errors

        for ((balance=0; balance<18; balance++))
        do

            #name of trained model
            trainModel="train_0${fold}_${balance}.model"
            echo "Testing model from fold ${fold} + balance ${balance} on element ${element}"

    	   # This tests on a pre-trained model
    	   java $CLASSIFIER -l $trainModel -T $testARFF -p 0 > test2_${fold}_${element}_${balance}.out
           #will give test2_0_0_0.out format output files

        done

        echo "Calculating TP/TN/FP/FN output for fold ${fold} for element ${element}"
        for file in ./test2_${fold}_${element}_*.out
        do
             $bindir/confusion_matrix_calc.py $file >> matrix_results_${fold}_${element}.txt
             #all the TP/FP/TN/FN results in a text file for that datapoint across all balancing runs

        done

        #get the mode from the matrix_results file to get the result for that element (TP/FP/TN/FN)
        $bindir/matrix_mode.py matrix_results_${fold}_${element}.txt >> mcc_data_${fold}.txt

        element=$(expr $element + 1)

    done

    #calculate MCC score for this fold
    echo "Calculating MCC score for fold ${fold}"
    $bindir/mcc_calc.py mcc_data_${fold}.txt >> mcc_scores.txt

    mv *.arff *.errors *${fold}.csv test_files_${fold}

done

#get mean of the 10 mcc scores
echo "The MCC mean across the 10 folds is equal to:"
$bindir/mcc_scores.py mcc_scores.txt 




