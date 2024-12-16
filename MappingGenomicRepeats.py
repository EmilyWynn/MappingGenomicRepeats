#!/usr/bin/env python3

#This python script will identify genomic regions of interest likely caused by the presence of repeats
#Before running this script, make a folder with the sample and library ID
#In Geneious, export the "Contig" file generated when remapping the falcon reads to the consensus
#Export this file to the new folder and save it as a ".bam" file
#Check the box to include the reference fasta when exporting


import os, sys, csv

#This section of code takes the bam and fasta files and indexes them
#Then runs them through mpileup, which compiles all SNPs and indels found in the alignment
#bcftools query changes the format into a more python friendly version

os.system("samtools faidx Contig.fasta")
os.system("samtools index -b Contig.bam")
os.system("samtools mpileup -v -f Contig.fasta Contig.bam | bcftools call -mv -O v -o pileup.vcf")
os.system("bcftools query -f '%POS %REF %ALT %INFO/IMF\n' pileup.vcf > temp.csv")

#This section of code takes the bcftools query output and analyzes it for regions of interest
#It will identify regions of many nearby SNPs and output the beginning of the region
#It will also identify regions where there is an indel in a high percentage of those reads
#Both of these are signs that there is a repeat nearby that is causing incorrect mapping
#By default it looks for SNPs that are at least 20bp away. This can be changed by changing the SNPwindow value
#By default it looks for INDELS that are at a frequency of 0.75. This can be changed by changing the INDELfreq value

file = open("regions.txt","w")
with open('temp.csv') as csvinput:
	reader = csv.reader(csvinput, delimiter=' ')
	pos = []
	ref = []
	alt = []
	freq = []
	SNPwindow = 20
	INDELfreq = 0.75
	for row in reader:
		pos.append(row[0])
		ref.append(row[1])
		alt.append(row[2])
		freq.append(row[3])
	snpregion = [0] * len(pos)
	for x in range(len(pos)-SNPwindow):
		if len(ref[x]) == 1 and len(alt[x]) == 1:
			for y in range(SNPwindow):
				if str(int(pos[x])) == str(int(pos[x+y])):
					snpregion[x] = 1
				else:
					continue
		else:
			snpregion[x] = '0'
	for x in range(len(snpregion)-3):
		if str(snpregion[x]) == '0' and str(snpregion[x+1]) == '1' and str(snpregion[x+2]) == '1' and str(snpregion[x+3]) == '1' and str(snpregion[x-1]) == '0':
			print(pos[x+1])
			file.write(pos[x+1])
			file.write("\n")
	for x in range(len(freq)):
		if freq[x] == '.':
			continue
		if float(freq[x]) > INDELfreq and len(alt[x])-len(ref[x]) > 2 or float(freq[x]) > INDELfreq and len(ref[x])-len(alt[x]) > 2:
			print(pos[x])
			file.write(pos[x])
			file.write("\n")
os.system("rm temp.csv")
