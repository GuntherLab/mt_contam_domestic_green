#!/usr/bin/env python

import operator
from optparse import OptionParser

parser = OptionParser()

parser.add_option('-i',action='store',type='string',dest='infile',help='Input alignment in FASTA format with reference as first sequence')
parser.add_option('-o',action='store',type='string',dest='outfile',help='output file name',default='MAFfile.txt')

(options, args) = parser.parse_args()
infile=options.infile
outfile=options.outfile


f=open(infile)
lines=f.readlines()
f.close()

seqs=[]
s=''

for l in lines[1:]:
	if l[0]=='>':
		seqs.append(s)
		s=''
		continue
	s+=l.rstrip()

nucs=['A','C','G','T','-']
pos=0
indices=[]
ind_pos={}

for i,x in enumerate(seqs[0]):
	if x!='-':
		pos+=1
		ind_pos[i]=pos
		indices.append(i)
		
#print(seqs[0][300:330])

f=open(outfile,'w')

for i in indices:
	counts={}
	for x in nucs:
		counts[x]=0



	for s in seqs[1:]:
		if s[i] in nucs:
			counts[s[i]]+=1

	cons_allele=max(counts.iteritems(),key=operator.itemgetter(1))[0]
	cons_af=float(max(counts.values()))/sum(counts.values())

	if 1: #cons_allele!='-':
		pos+=1
		#if cons_af>=0.99:
		#	private=1
		#else:
		#	private=0
		f.write('%s\t%s\t%s\n' %(ind_pos[i],cons_allele,cons_af))

f.close()
