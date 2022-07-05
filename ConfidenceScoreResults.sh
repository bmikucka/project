#!/bin/sh


# Scripts
bindir=/serv/data/people/bmikucka/bin

#for each fold
for ((fold=0; fold<10; fold++))
do
	#for each of the individual datapoint files for that fold - for each file test_0${fold}_${element}.csv
    element=0
    for file in ./test_0${fold}_*.csv #this is for one element
    do
		#for each balancing run result of that element- get confidence score
		#!!! CORRECT LOCATION OF OUT FILES
		for file in ./test2_${fold}_${element}_*.out #this is for each balance run for one element
        do
             $bindir/confidence_score.py $file >> confidence_score_${fold}_${element}.txt
        done
        
		#average confidence scores to get TP/FP/TN/FN predicition (mcc_data_0.txt a list of FP/TP... to use for MCC calculation)
		$bindir/datapoint_prediction.py confidence_score_${fold}_${element}.txt >> mcc_data_${fold}_cs.txt

		element=$(expr $element + 1)
	done

	#get mcc score for that fold
	echo ""
	echo "Calculating MCC score for fold ${fold}"
	$bindir/mcc_calc.py mcc_data_${fold}_cs.txt >> mcc_scores_cs.txt

done

#get overall mcc mean
echo ""
echo ""
echo ""
echo "The MCC mean across the 10 folds is equal to:"
$bindir/mcc_scores.py mcc_scores_cs.txt 