#!/usr/bin/env python3

'''
Author: Dario Galanti, April 2020
AIM: Filtering and Imputation of NAs in union_bedGraph files. This script calculates a beta distribution from each marker and uses it to randomly sample NAs at that position.
AIM: Such imputation method is minimizing the effect of NAs on the position, it is suitable to estimate few NAs per position without removing all positions with at least 1 NA
RUN: sbatch --partition test --cpus-per-task 4 --mem 24G --time 06:00:00 --wrap "python3 beta_impute_ewas.py CpG.unfilter.bed"
Note: Positions with all 1s or all 0s are imputed with a 1 or a 0 respectively. Beta distribution is not suitable to model these positions.
NB: Set filter_NA in line 19 and filter_SD in line 20!! If willing to keep also lines with SD==0 a different script should be used, ask Dario

Dependencies: numpy is required

'''

import sys
from statistics import mean 
from statistics import pstdev 
import numpy as np
import argparse

argparser = argparse.ArgumentParser()
argparser.add_argument('--filter_NA', type=double)
argparser.add_argument('--filter_SD', type=double)
ARGS = argparser.parse_args()

#print(ARGS.text)

fin = open(sys.argv[1], "r")
fout = open("beta_imputed_" + sys.argv[1], "w")

#filter_NA=${params.filter_NA}									#Max percentage of NAs of each line to be processed. Lines with more NAs are excluded.
#filter_SD=${params.filter_SD}								#Min standard dev of each line to be processed. Cannot be zero!!!

## 1: Loop through lines with NAs < filter_NA.
line_num = 0
proc_lines = 0
for line in fin:
	line = line.strip()
	splitline = line.split("\t")
	line_num += 1
	pos = splitline[:3]
	val = splitline[3:]
	spls = len(val)
	max_NA = spls * filter_NA						#Calculate max number of NAs
	
## 2: Calculate position average and st.dev
	if line_num > 1 and line.count('NA') < max_NA:
		v = [i for i in val if str(i) != 'NA']
		v = [float(i) for i in v]
		av = mean(v)
		sd = pstdev(v)
		
		if sd > filter_SD:
			proc_lines += 1
## 3: Fit the beta distribution (manual calculation of a and b parameters)
			if 0.001 < av < 0.999:					#Exclude positions with all zeros or all ones and sd = 0 (sd = 0 won't allow a and b calculation)
				a = av**2 * ((1 - av) / sd**2 - (1 / av))
				b = a * (1 / av - 1)
## 4: Randomly sample values from the beta distribution with "numpy" and substitute them to NAs
				for x in range(spls):
					if val[x] == "NA":
						val[x]=	np.random.beta(a,b)
						val[x]=	'{:.2f}'.format(val[x])
					
			elif av >= 0.999:						#Impute positions with all "1s"
				for x in range(spls):
					if val[x] == "NA":
						val[x]= '{:.2f}'.format(1.00)
			else:									#Impute positions with all "0s"
				for x in range(spls):
					if val[x] == "NA":
						val[x]= '{:.2f}'.format(0.00)
			pos.extend(val)							#Append val list to pos list
			print(*pos, sep="\t", file=fout)
	elif line_num < 2:
		print(line, sep="\t", file=fout)		 	#Print headers

fin.close()
fout.close()

#print(str(sys.argv[1]) + " contains " + str(spls) + " individulas")
#print("Out of " + str(line_num) + " initial markers, " + str(proc_lines) + " were kept after filtering and imputation")
