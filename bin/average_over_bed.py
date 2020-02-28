#!/usr/bin/env python

'''
Title: average_over_bed.py
Date: 20191012
Author: Adam Nunn
Description:
	This script iterates over a positions file generated by bedtools unionbedg
	and a corresponding regions file, to give the average values of all positions
	contained within each region. Only works with position sorted bed files with
	distinct, non-overlapping regions, such as the DMRs from Metilene.

List of functions:
	main()
	moveToNextRegion()

Procedure:
	1. Open the regions and positions files for reading
	2. Iterate through the positions file
	3. On each position, test if we need to move to the next region
	4. On each position, evaluate or skip to next

Usage:
	average_over_bed.py [REGIONS] [POSITIONS]
eg. average_over_bed.py regions.bed positions.bed
'''

###################
## INIT ENVIRONMENT

import sys


#################
## BEGIN __MAIN__
def main(REGIONS,POSITIONS):

	# 1) Open the regions and positions files for reading
	with open(REGIONS,'r') as regions, open(POSITIONS,'r') as positions:

		# setup vars for initial region
		try: region = next(regions)
		except StopIteration: 
			print("No regions detected in bedfile")
			raise SystemExit(0)
		
		region = region.rstrip()
		region = region.split("\t")

		Chr = region[0]
		Start = int(region[1])
		End = int(region[2]) + 1

		# define booleans for header and region detection
		Head = True
		Found = False

		# 2) Iterate through the positions file
		for position in positions:

			# split the position line
			position = position.rstrip()
			position = position.split("\t")

			# stage variables and print the header
			if Head:

				# establish initial variables
				length = len(position) - 2
				totals = [(0, 0) for i in range(0,length)]
				#count = 0

				# print header and skip to next position
				position.insert(1, "start")
				print("\t".join(position))
				Head = False
				continue

			##############################################
			# establish current position
			current_Chr = position[0]
			current_Pos = int(position[1])

			# 3) On each position, test if we need to move to the next region
			if Found and ((current_Chr != Chr) or ((current_Chr == Chr) and (current_Pos > End+1))):

				# print averages for current region
				#averages = [str(format(n/count, '.5f')) for n in totals]
				averages = [str(format(n[0]/n[1], '.5f')) if n[1] != 0 else "NA" for n in totals]
				print("{}\t{}\t{}\t{}".format(Chr, Start, End-1, "\t".join(averages)))

				# move the region ONCE
				region, Chr, Start, End = moveToNextRegion(regions, Chr, Start, End)
				if region is "kill": break

				# reset counts
				totals = [(0, 0) for i in range(0,length)]
				Found = False

			# 4) On each position, evaluate current position or skip to next
			if (current_Chr == Chr) and (current_Pos in range(Start, End+1)):

				totals = [(totals[i][0],totals[i][1]) if position[i+2] == "NA" else (totals[i][0]+float(position[i+2]),totals[i][1]+1) for i in range(0, length)]
				Found = True

			# we already moved the region, so now move position
			else: continue


## END OF __MAIN__
##################


###################
## DEFINE FUNCTIONS

# move to the next region until the next region contains position
def moveToNextRegion(regions, Chr, Start, End):

	# move the region
	try: region = next(regions)
	except StopIteration: return "kill", Chr, Start, End

	region = region.rstrip()
	region = region.split("\t")
	Chr = region[0]
	Start = int(region[1])
	End = int(region[2]) + 1
	
	return region, Chr, Start, End


## END OF FUNCTIONS
###################

#############
## RUN SCRIPT

# run main()
if __name__ == '__main__':
	main(sys.argv[1],sys.argv[2])

## END OF SCRIPT
################