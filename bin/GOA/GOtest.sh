#!/bin/bash

hostname
date

infile=$1
species=$2
OverTree=$3
Over=$4
UnderTree=$5
Under=$6
GOcat=$7
Pval=$8

R CMD BATCH --vanilla "--args $infile $species $OverTree $Over $UnderTree $Under $GOcat $Pval" /home/cansu/Desktop/GO_Wordle/GOtest_bin.R 

if [ -f $Over ]
then
	python /home/cansu/Desktop/GO_Wordle/createWordle.py $Over "Enhancement"
else
	echo "No data for $Over"
fi

if [ -f $Under ]
then
	python /home/cansu/Desktop/GO_Wordle/createWordle.py $Under "Suppression"
else
	echo "No data for $Under"
fi

date
