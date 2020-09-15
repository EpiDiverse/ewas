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

<<<<<<< HEAD
R CMD BATCH --vanilla "--args $infile $species $OverTree $Over $UnderTree $Under $GOcat $Pval" ${baseDir}bin/GOA/GOtest_bin.R 

if [ -f $Over ]
then
	python ${baseDir}bin/GOA/createWordle.py $Over "Enhancement"
=======
R CMD BATCH --vanilla "--args $infile $species $OverTree $Over $UnderTree $Under $GOcat $Pval" ${baseDir}/bin/GOA/GOtest_bin.R 

if [ -f $Over ]
then
	python ${baseDir}/bin/GOA/createWordle.py $Over "Enhancement"
>>>>>>> 3ed2bddd5686fcd937992d9f55274d29ea7fd703
else
	echo "No data for $Over"
fi

if [ -f $Under ]
then
<<<<<<< HEAD
	python ${baseDir}bin/GOA/createWordle.py $Under "Suppression"
=======
	python ${baseDir}/bin/GOA/createWordle.py $Under "Suppression"
>>>>>>> 3ed2bddd5686fcd937992d9f55274d29ea7fd703
else
	echo "No data for $Under"
fi

date
