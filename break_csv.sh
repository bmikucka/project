#!/bin/sh

file=$1
fold=$2

#get number of elements
num=$(cat $file | wc -l)
num=$(expr $num - 1) #without the headline

for ((element=0; element<$num; element++))
do
	head -1 file > test_0${fold}_${element}.csv
	line_num=$(expr $element + 2) #starts from 0 and avoid head line
	sed '${line_num}!d' file >> test_0${fold}_${element}.csv

done
