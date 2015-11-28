import re
import numpy as np
import json
import csv
from collections import Counter


"""Setting counter"""
#Output: empty json structure for codon counter
def createCodonCounter():
	template={
		'F':{'TTT':0,'TTC':0},
		'L':{'TTA':0,'TTG':0,'CTT':0,'CTC':0,'CTA':0,'CTG':0},
		'I':{'ATT':0,'ATC':0,'ATA':0,},
		'M':{'ATG':0},
		'V':{'GTT':0,'GTC':0,'GTA':0,'GTG':0},
		'S':{'TCT':0,'TCC':0,'TCA':0,'TCG':0,'AGT':0,'AGC':0},
		'P':{'CCT':0,'CCC':0,'CCA':0,'CCG':0},
		'T':{'ACT':0,'ACC':0,'ACA':0,'ACG':0},
		'A':{'GCT':0,'GCC':0,'GCA':0,'GCG':0},
		'Y':{'TAT':0,'TAC':0},
		'*':{'TAA':0,'TAG':0,'TGA':0},
		'H':{'CAT':0,'CAC':0},
		'Q':{'CAA':0,'CAG':0},
		'N':{'AAT':0,'AAC':0},
		'K':{'AAA':0,'AAG':0},
		'D':{'GAT':0,'GAC':0},
		'E':{'GAA':0,'GAG':0},
		'C':{'TGT':0,'TGC':0},
		'W':{'TGG':0},
		'R':{'CGT':0,'CGC':0,'CGA':0,'CGG':0,'AGA':0,'AGG':0},
		'G':{'GGT':0,'GGC':0,'GGA':0,'GGG':0}
	}
	json.dumps(template)
	return template


"""Opening files"""

with open("orf_trans.fasta protein","rb")as f:
#with open("orf_trans.fasta protein","rb")as f:
	my_list=[line.rstrip('\n') for line in f]

with open("orf_coding.fasta_DNA","r") as f2:
	my_list2=[line.rstrip('\n') for line in f2]

"""Manipulating files"""
#Input: file
#Output: list of protein name, list of protein sequence
def separate(MyFile):
	var=0
	str1=None
	string=None
	temp1=[]
	temp2=[]

	while var < len(MyFile):
		plusOne=var+1
		#print plusOne
		if MyFile[var][0] is (">"):
			splitted=MyFile[var].split()[:3]
			str1=" ".join(splitted)
			temp1.append(str1)
		else:
			if string is None:
				string=MyFile[var]
			else:
				string=string+MyFile[var]

			if plusOne <(len(MyFile)):
				if MyFile[plusOne][0] is (">"):
					temp2.append(string)
					string=None

			if plusOne == (len(MyFile)):
				temp2.append(string)
				string=None
				return temp1,temp2

		var=var+1
	return temp1,temp2

#Input: name of sequences, sequences
#Output: list of name at 0 position and following the infomation
def merge(name,sequence):
	my_list=[]
	for i in range(len(name)):
		my_list.append([name[i],sequence[i]])
	return my_list


"""Find homopeptides & homocodons"""

#Input: a sequence
#Output: list of homopeptides [length, amino acid that's repeating,starting position]
def isHomopeptide(sequence):
	temp=[]
	pattern=re.compile(r'((\w)\2{2,})')
	for repeat in pattern.finditer(sequence):
		start=repeat.start()
		group=repeat.group()
		aa=group[0]
		length=len(group)
		temp.append([length,aa,start])
	return temp

#Input: list of all sequence in file
#Output: list of all homopeptides in all sequences
def homopeptide(sequence_list):
	result=[]

	for i in range(len(sequence_list)):
		result.append(isHomopeptide(sequence_list[i]))

	return result

def isHomocodon(dnaSeq,length):
	initial=dnaSeq[0]

	for i in range(length):
		codon=dnaSeq[i]
		if codon != initial:
			return False

	return True

#Inputs: list of the homopeptides(with details), DNA sequence, codon counter
#Outputs: codon counter
def findHomocodons(protein_sequence,DNA_sequence):
	codons=[]
	same=[]
	notSame=[]
	CodonCounter=createCodonCounter()
	homocodonCounter=createCodonCounter()
	noHMCCOunter=createCodonCounter()

	for i in range(len(protein_sequence)):
		# print protein_sequence[i][0]
		for j in range(len(protein_sequence[i])):
			# print ("protein_sequence[i][j]")
			# print protein_sequence[i][j]
			aminoA=protein_sequence[i][j][1]
			length=protein_sequence[i][j][0]
			
			startPosition=protein_sequence[i][j][2]
			start=(startPosition* 3)
			
			endPosition=startPosition+length
			end=(endPosition*3)
			
			#store=codon[0:3]
			#check=True

			codon=[]

			for x in range(length):
				codon.append(DNA_sequence[i][start+(x*3):start+((x+1)*3)])

			check=isHomocodon(codon,length)

			if check:
				for y in range(length):
					CodonCounter[aminoA][codon[y]]=CodonCounter[aminoA][codon[y]]+1
					homocodonCounter[aminoA][codon[y]]=homocodonCounter[aminoA][codon[y]]+1
				same.append([i,length,aminoA,startPosition,codon[0]])
			else:
				for y in range(length):
					CodonCounter[aminoA][codon[y]]=CodonCounter[aminoA][codon[y]]+1
					noHMCCOunter[aminoA][codon[y]]=noHMCCOunter[aminoA][codon[y]]+1
				notSame.append([i,length,aminoA,startPosition,codon])

	return CodonCounter,homocodonCounter,noHMCCOunter,same,notSame

def findNoHomopeptides(homopeptides,DNASequence,proteinSequence):
	counter1=createCodonCounter()  # homopeptides
	counter2=createCodonCounter()  # others


	for i in range(len(homopeptides)):
		position=[]
		protein=proteinSequence[i]

		for j in range(len(homopeptides[i])):
			# print (homopeptides[i][j])
			for k in range(homopeptides[i][j][0]):
				position.append(homopeptides[i][j][2]+k)
		for x in range(len(protein)):
			aminoacid=protein[x]
			DNA=DNASequence[i][x*3:(x+1)*3]

			if x in position:
				counter1[aminoacid][DNA] += 1
			else:
				counter2[aminoacid][DNA] += 1
	
	return counter1,counter2


"""Find gaps"""

def findGaps(homopeptides,proteinSeq,aa):
	sequenceNumber=len(homopeptides)
	mylist=[]
	gapInfo=[]


	for i in range(sequenceNumber):
		j=0

		while (j+1)<len(homopeptides[i]):
			aminoA=homopeptides[i][j][1]
			length=homopeptides[i][j][0]
			startPosition=homopeptides[i][j][2]

			aminoA2=homopeptides[i][j+1][1]
			length2=homopeptides[i][j+1][0]
			startPosition2=homopeptides[i][j+1][2]
			endPosition2=startPosition2+length2

			if aminoA == aminoA2:
				if aminoA == aa:
					if startPosition2 < startPosition+length+3:
						gap=proteinSeq[i][startPosition+length:startPosition2]
						length3=len(gap)
						totalLength=length+length2+length3
						gapProportion=length3/float(totalLength)
						mylist.append([i,homopeptides[i][j],gap,homopeptides[i][j+1]])
						gapInfo.append([i,length3,totalLength,gapProportion])
		
			j=j+1		

	return mylist,gapInfo

def Gapcounter(gaplist,DNASequence):
	counter=createCodonCounter()

	for i in range(len(gaplist)):
		index=gaplist[i][0]
		info=gaplist[i][1]
		amino=gaplist[i][2]
		position=info[2]+info[0]
		if len(amino)==1:
			DNA=DNASequence[index][position*3:(position+1)*3]
			counter[amino][DNA]+=1
		else:
			a=amino[0]
			aa=amino[1]
			DNA=DNASequence[index][position*3:(position+1)*3]
			counter[a][DNA]+=1
			position2=position+1
			DNA2=DNASequence[index][position2*3:(position2+1)*3]
			counter[aa][DNA2]+=1
	return counter

"""Consecutive homopeptides"""

def consecutive(homopeptides):
	var=len(homopeptides)
	my_list=[]

	for i in range(var):
		for j in range(len(homopeptides[i])-1):
			length=homopeptides[i][j][0]
			startPosition=homopeptides[i][j][2]
			value=startPosition+length
			# print(length,startPosition,value)

			startPosition2=homopeptides[i][j+1][2]
			if value==startPosition2:
				my_list.append([i,homopeptides[i][j],homopeptides[i][j+1]])
	return my_list

def consecutiveSeq(my_list,dnaSeq):
	f=open('result.txt','w')

	for i in range(len(my_list)):
		entry=my_list[i]
		startPosition=entry[1][2]
		endPosition=entry[2][2]+entry[2][0]
		geneIndex=entry[0]
		DNASequence=dnaSeq[geneIndex][startPosition*3:endPosition*3]
		print(my_list[i],DNASequence)

	return

def countPairs(my_list):
	Ndict={}
	Qdict={}

	for i in my_list:
		first=i[1]
		second=i[2]
		a1=first[1]
		a2=second[1]

		if a1 == 'N':
			pair='N'+a2
			if pair in Ndict:
				Ndict[pair]+=1
			else:
				Ndict[pair]=1
		elif a2 == 'N':
			pair='N'+a1
			if pair in Ndict:
				Ndict[pair]+=1
			else:
				Ndict[pair]=1
		elif a1 == 'Q':
			pair=a1+a2
			if pair in Qdict:
				Qdict[pair]+=1
			else:
				Qdict[pair]=1
		elif a2 == 'Q':
			pair=a2+a1
			if pair in Qdict:
				Qdict[pair]+=1
			else:
				Qdict[pair]=1
	with open("Npair.csv","wb") as csvfile:
		f=csv.writer(csvfile)
		f.writerow(['pair','count'])
		for pair in Ndict:
			count=Ndict[pair]
			f.writerow([pair,count])
	with open("Qpair.csv","wb") as csvfile:
		f=csv.writer(csvfile)
		f.writerow(['pair','count'])
		for pair in Qdict:
			count=Qdict[pair]
			f.writerow([pair,count])

	return Ndict,Qdict


"""Analysis"""

#Input: the counter for codons
#Output: the proportion for each codon in every amino acids 
def proportion(template,filename,name):
	with open(filename,'wb') as csvfile:
		f=csv.writer(csvfile)
		f.writerow(['amino acid','codon',name])
		for AminoAcid in template:
			total = 0
			#print AminoAcid
			for codons in template[AminoAcid]:
				total= total+template[AminoAcid][codons]

			for codons in template[AminoAcid]:
				if total is not 0:
					value=float(template[AminoAcid][codons])/float(total)
					#print codons,total,value
					f.writerow([AminoAcid, codons, value])
				else:
					f.writerow([AminoAcid, codons, total])
	return

#Input: the counter for codons
#Output: the maximum codon for each amino acids
def findMax(template,filename):
 	maxValue=0
 	maxCodon=None
 	result=[]
 	for AminoAcid in template:
 		for codons in template[AminoAcid]:
	 		if template[AminoAcid][codons]> maxValue:
	 			maxValue=template[AminoAcid][codons]
	 			temp=template[AminoAcid]
	 			proportion=float(maxValue)/sum(temp.values())
	 			maxCodon=codons
	 	result.append([AminoAcid,maxCodon,proportion])
	 	maxValue=0
	 	maxCodon=None

 	with open(filename,'wb') as csvfile:
 		f=csv.writer(csvfile)
		f.writerow(['amino acid','codon','proportion'])
		for i in range(len(result)):
			aa=result[i][0]
			codon=result[i][1]
			mc=result[i][2]
			f.writerow([aa,codon,mc])
	return

def toCSV(template,filename):
	with open(filename,'wb') as csvfile:
		f=csv.writer(csvfile)
		f.writerow(['amino acid','codon'])
		for aa in template:
			total=0
			for codon in template[aa]:
				k=template[aa][codon]
				total+=k

			f.writerow([aa,total])


proteinName,proteinSequence=separate(my_list)
proteinHomopeptides=homopeptide(proteinSequence)
# TotalHomopeptides=merge(proteinName,proteinHomopeptides)

DNAName,DNASequence=separate(my_list2)
# counter,Hcounter,NHcounter,same,notsame=findHomocodons(proteinHomopeptides,DNASequence)
# homocounter,othercounter=findNoHomopeptides(proteinHomopeptides,DNASequence,proteinSequence)
# proportion(homocounter,"HOMOPEPTIDES.csv","homopeptides")
# proportion(othercounter,"OTHERS.csv","others")

# findMax(homocounter,"homopeptidesMC.csv")
# findMax(othercounter,"otherMC.csv")

# aminoAcid='Q'
# aminoAcid='N'
# gaplist,gapinfo=findGaps(proteinHomopeptides,proteinSequence,aminoAcid)
# for i in range(len(gaplist)):
#  	print(gaplist[i])
#  	# print(gapinfo[i])
# print (len(gaplist))

# gapC=Gapcounter(gaplist,DNASequence)
# # toCSV(gapC,"Qgapcount.csv")
# toCSV(gapC,"Ngapcount.csv")
# with open('HOMOPEPTIDES.csv','wb') as csvfile:
# 	f=csv.writer(csvfile)
# 	f.writerow(['amino acid','codon','homopeptides'])
# 	for i in homocounter:
# 		for j in homocounter[i]:
# 			k=homocounter[i][j]
# 			f.writerow([i,j,k])

consec=consecutive(proteinHomopeptides)
# for i in range(len(consec)):
#  	print(consec[i])
#  	#print(gapinfo[i])

# print(len(consec))

#consecutiveSeq(consec,DNASequence)

countPairs(consec)

# with open('OTHERS.csv','wb') as csvfile:
# 	f=csv.writer(csvfile)
# 	f.writerow(['amino acid','codon','others'])
# 	for i in othercounter:
# 		for j in othercounter[i]:
# 			k=othercounter[i][j]
# 			f.writerow([i,j,k])
# propensities = {'A' :-0.396490246, 'C' : 0.415164505, 'D' : -1.276997939, 'E' : -0.605023827, 'F' : 0.838732498,
#                'G' : -0.039220713, 'H': -0.278573356, 'I' : 0.813697862, 'K': -1.576748587, 'L' : -0.040005335,
#                'M' : 0.673729095, 'N' : 0.080295334, 'P' : -1.197447496, 'Q' : 0.069168387, 'R' : -0.405858577,
#                'S' : 0.133912418, 'T' : -0.11457038, 'V' : 0.813697862, 'W' : 0.666735081, 'Y' : 0.77865336}
# with open('PAPApropensities.csv','wb') as csvfile:
# 	f=csv.writer(csvfile)
# 	f.writerow(['amino acid','propensities'])
# 	for aa in propensities:
# 		prop=propensities[aa]
# 		f.writerow([aa,prop])
