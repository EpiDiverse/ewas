#!/usr/bin/env python
from __future__ import division
import sys
import numpy
fin=open(sys.argv[1],'r')
f_tend=sys.argv[2]
if f_tend=="Enhancement":
    fout=open(sys.argv[1].split('.')[0]+'_Enhancement_forWordleNew.txt','w')
if f_tend=="Suppression":
    fout=open(sys.argv[1].split('.')[0]+'_Suppression_forWordleNew.txt','w')
#fsup=open(sys.argv[1].split('.')[0]+'_suppression.txt','w')
#fenh=open(sys.argv[1].split('.')[0]+'_enhancement.txt','w')
line=fin.readline()
lines=fin.readlines()
lines=[x for x in lines if x.strip().split('\t')[1:]!=['NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA']]
newLine=[]
supLine=[]
enhLine=[]
for line in lines:
	p=float(line.split('\t')[8][0:-1])
	p=abs(numpy.log10(p))
	if f_tend=="Enhancement":		
		if p<4:
			newLine.append(line.split('\t')[7][1:-1]+":"+str(p)+':A2BC13\n') #Kermit
		if p>=4:
			newLine.append(line.split('\t')[7][1:-1]+":"+str(p)+':385E0F\n') #terreverte
#		enhLine.append(line.split('\t')[0]+'\n')
	if f_tend=="Suppression":
		if p<4:
			newLine.append(line.split('\t')[7][1:-1]+":"+str(p)+':FF8C69\n') #salmon1
		if p>=4:
			newLine.append(line.split('\t')[7][1:-1]+":"+str(p)+':8C1717\n') #SCARLET
#		supLine.append(line.split('\t')[0]+'\n')
	#newLine.append(line.split('\t')[0]+'\n')
fout.write(''.join(newLine))
#fenh.write(''.join(enhLine))
#fsup.write(''.join(supLine))
