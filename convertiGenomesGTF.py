#!/usr/bin/env python
# Author Fong Chun Chan <fongchunchan@gmail.com>
import argparse

parser = argparse.ArgumentParser(description='This script is used to converting the iGenomes GTF file from Ensembl so that we have a chr added (for chrMT it is converted to M) to and we remove the e chromosomes')
parser.add_argument('inGTF', action='store', help='The iGenomesGTF File downloaded from either Tophat or Cufflinks website')
args = parser.parse_args()

inFh = open(args.inGTF, 'rb')
chr_list = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','21','X','Y','MT']
for line in inFh:
	vals = line.rstrip().split("\t")	
	if vals[0] in chr_list:
		if vals[0] == 'MT':
			print "chrM" + "\t" + "\t".join(vals[1:len(vals)])
		else :
			print "chr" + vals[0] + "\t" + "\t".join(vals[1:len(vals)])

	
