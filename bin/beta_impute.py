#!/usr/bin/env python3

'''
Author: Dario Galanti, April 2020 (modified by Nilay Can)
AIM: Filtering and Imputation of NAs in union_bedGraph files. This script calculates a beta distribution from each marker and uses it to randomly sample NAs at that position.
AIM: Such imputation method is minimizing the effect of NAs on the position, it is suitable to estimate few NAs per position without removing all positions with at least 1 NA
RUN: sbatch --partition test --cpus-per-task 4 --mem 24G --time 06:00:00 --wrap "python3 beta_impute_ewas.py CpG.unfilter.bed"
Note: Positions with all 1s or all 0s are imputed with a 1 or a 0 respectively. Beta distribution is not suitable to model these positions.
NB: Set filter_NA in line 19 and filter_SD in line 20!! If willing to keep also lines with SD==0 a different script should be used, ask Dario

Dependencies: numpy is required

'''
#############
## INIT ENVIRONMENT

import sys
from statistics import mean 
from statistics import pstdev 
import numpy as np
import argparse

'''
argparser = argparse.ArgumentParser()
argparser.add_argument('--filter_NA', type=double)
argparser.add_argument('--filter_SD', type=double)
ARGS = argparser.parse_args()
'''
#print(ARGS.text)

#############
#DEFINE __MAIN__
def main(fin,fout):
    
    #1) Open the unfiltered file for reading and writing
    with open(sys.argv[1], "r") as fin, open("beta_imputed_" + sys.argv[1], "w") as fout:
        
        
        ## 2: Loop through lines with NAs < filter_NA.
        line_num = 0
        proc_lines = 0
        for line in fin:
            line = line.strip()
            splitline = line.split("\t")
            line_num += 1
            pos = splitline[:3]
            val = splitline[3:]
            spls = len(val)
            max_NA = spls * filter_NA                       #Calculate max number of NAs
            
            ## 3: Calculate position average and st.dev
            if line_num > 1 and line.count('NA') < max_NA:
                v = [i for i in val if str(i) != 'NA']
                v = [float(i) for i in v]
                av = mean(v)
                sd = pstdev(v)
                
                if sd > filter_SD:
                    proc_lines += 1
                    
                    
                    ## 4: Fit the beta distribution (manual calculation of a and b parameters)
                    if 0.001 < av < 0.999:                  #Exclude positions with all zeros or all ones and sd = 0 (sd = 0 won't allow a and b calculation)
                        a = av**2 * ((1 - av) / sd**2 - (1 / av))
                        b = a * (1 / av - 1)
                        
                        ## 5: Randomly sample values from the beta distribution with "numpy" and substitute them to NAs
                        for x in range(spls):
                            if val[x] == "NA":
                                val[x]= np.random.beta(a,b)
                                val[x]= '{:.2f}'.format(val[x])
                                
                                
                    elif av >= 0.999:                       #Impute positions with all "
                        for x in range(spls):
                            if val[x] == "NA":
                                val[x]= '{:.2f}'.format(1.00)
                    else:                                   #Impute positions with all "0s"
                        for x in range(spls):
                            if val[x] == "NA":
                                val[x]= '{:.2f}'.format(0.00)
                    pos.extend(val)                         #Append val list to pos list
                    print(*pos, sep="\t", file=fout)
                elif line_num < 2:

## END OF __MAIN__
##################

#############
## RUN SCRIPT

# define argparse
#usage = 'read an ulfiltered file and impute missing values with beta distribution.'

        parser = argparse.ArgumentParser()
	
	parser.add_argument('infile', metavar='in.txt', 'the path to the unfiltered file')
	parser.add_argument('outfile', metavar='out.bed', 'the path to the beta imputed file')
	args = parser.parse_args()

# call main()
if __name__ == '__main__':
	main(args.infile,args.outfile)

## END OF SCRIPT
################

#print(line, sep="\t", file=fout)
#fin.close()
#fout.close()

#print(str(sys.argv[1]) + " contains " + str(spls) + " individulas")
#print("Out of " + str(line_num) + " initial markers, " + str(proc_lines) + " were kept after filtering and imputation")
