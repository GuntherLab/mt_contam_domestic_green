#!/usr/bin/env python
import sys
import operator
import math
from optparse import OptionParser


#samtools mpileup -B -q 30 -Q 30 -r MT -f $ref_seq $bamname | ./mt_contamination.py -s $MAF_file -f 0.05
#output: point estimate, informative sites, consensus alleles, total alleles, lower end CI, upper end CI


parser = OptionParser()

parser.add_option('-s',action='store',type='string',dest='fname',help='Sites MAF file')
parser.add_option('-f',action='store',type='float',dest='MAF_private',help='MAF below which sites are considered nearly private (default: 0.05)',default=0.05)

(options, args) = parser.parse_args()
fname=options.fname
MAF_private=options.MAF_private

f=open(fname)
lines=f.readlines()
f.close()

global_consensus={}
global_MAF={}

for l in lines:
	split=l.rstrip().split('\t')
	global_consensus[int(split[0])]=split[1]
	global_MAF[int(split[0])]=float(split[2])

nucs=['A','G','C','T','*']

consensus_all=0.0
all_all=0.0
sites=0
site_pos=[]

for l in sys.stdin:
	split=l.rstrip().split('\t')
	if int(split[3])>=10 and (1.0-global_MAF[int(split[1])])<=MAF_private:

		#if '-' in split[4] or '+' in split[4]:
		#	continue

		ref=split[2]
		ref_count=split[4].count(',')+split[4].count('.')
		counts={}

		for n in nucs:
			counts[n]=split[4].upper().count(n)

		counts[ref]=ref_count
		total=sum(counts.values())
		cons_allele=max(counts.iteritems(),key=operator.itemgetter(1))[0]

		if cons_allele=='C' or cons_allele=='G':
			continue

		if cons_allele!=global_consensus[int(split[1])]:
			#print global_consensus[int(split[1])],l
			print(' '.join(str(i) for i in ['Site: ',total, counts.values(),split[1]]))
			consensus_all+=max(counts.values())
			all_all+=total
			sites+=1
			site_pos.append(int(split[1]))

if all_all>0:
	c=1-consensus_all/all_all
	if consensus_all==all_all:
		conf_low=0.0
		conf_up=(1-math.pow(0.05,1.0/all_all))*100
	else:
		conf_up=min(1,c+1.96*math.sqrt((c*(1-c))/all_all))*100.0
		conf_low=max(0,c-1.96*math.sqrt((c*(1-c))/all_all))*100.0
		
	print(' '.join(str(i) for i in ['Point estimate:',c*100,'informative sites:',sites,'consensus alleles:',int(consensus_all),'total alleles:',int(all_all),'lower end CI:',conf_low,'upper end CI:',conf_up,'Positions:',site_pos]))

				
