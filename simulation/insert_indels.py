#!/usr/bin/python

import sys
import random

START=0
END=1

def usage():
	print "%s fasta_reference_filename min_size max_size number_indels output_fasta_file output_log_file prob_of_insert[0.5]" % sys.argv[0]

def overlaps(l,c_interval):
	#check to see if any overlap	
	for interval in l:
		#thanks to gfb
		#if max(interval["start"],c_interval["start"])>=min(interval["end"],c_interval["end"]):
		if interval[START]<=c_interval[START] and c_interval[START]<=interval[END]:
			#overlaps
			return True
		elif interval[START]<=c_interval[END] and c_interval[END]<=interval[END]:
			#overlaps
			return True
		elif c_interval[START]<=interval[START] and interval[END]<=c_interval[END]:
			#overlaps
			return True
		elif interval[START]<=c_interval[START] and c_interval[END]<=interval[END]:
			#overlaps
			return True 
	return False
	


letters=['A','C','G','T']

def random_string(l):
	s=""
	r=random.Random()
	while len(s)!=l:
		s+=r.choice(letters)
	return s

if __name__=='__main__':
	if len(sys.argv)!=7 and len(sys.argv)!=8:
		usage()
		sys.exit(1)

	prob_of_insert=0.5
	if len(sys.argv)==8:
		prob_of_insert=float(sys.argv[7])
	
	in_filename=sys.argv[1]
	min_size=int(sys.argv[2])
	max_size=int(sys.argv[3])
	num_indels=int(sys.argv[4])
	out_filename=sys.argv[5]
	log_filename=sys.argv[6]


	#read in the reference
	in_file=open(in_filename,'rU')
	ref_name=in_file.readline().strip()
	assert(ref_name[0]=='>')
	lines=in_file.readlines()
	ref_seq="".join(map(lambda x : x.strip(), lines))	
	in_file.close()

	#think of some insertions
	indel_intervals=[]
	original_intervals=[]
	r=random.Random()
	j=0
	while j<num_indels:
		size=r.randrange(min_size,max_size+1)
		start=r.randrange(len(ref_seq)-size)
		candidate_interval=[start,start+size-1]
		if not overlaps(indel_intervals,candidate_interval) and ref_seq[start:start+size].upper()==ref_seq[start:start+size]:
			indel_intervals.append(candidate_interval[:])
			original_intervals.append(candidate_interval[:])
			j=j+1
	indel_intervals.sort()
	original_intervals.sort()

	previous=0
	ref_seq_windels=""
	for x in range(len(indel_intervals)):
		insert=r.random()<=prob_of_insert
		interval=indel_intervals[x]
		original_interval=original_intervals[x]
		if insert:
			ref_seq_windels=ref_seq_windels+ref_seq[previous:original_interval[START]]+random_string(original_interval[END]-original_interval[START]+1).lower()
			for y in range(x+1,len(indel_intervals)):
				indel_intervals[y][START]+=(interval[END]-interval[START]+1)
				indel_intervals[y][END]+=(interval[END]-interval[START]+1)
			interval.append("Insert")
			original_interval.append("Insert")
			previous=original_interval[START]
		else:
			ref_seq_windels=ref_seq_windels+ref_seq[previous:original_interval[START]]
			for y in range(x+1,len(indel_intervals)):
				indel_intervals[y][START]-=(interval[END]-interval[START]+1)
				indel_intervals[y][END]-=(interval[END]-interval[START]+1)
			interval.append("Deletion")
			original_interval.append("Deletion")
			previous=original_interval[END]+1
	ref_seq_windels=ref_seq_windels+ref_seq[previous:]

	#print out the log
	log_file=open(log_filename,'w')
	print >> log_file, "\n".join(map( lambda (x) : "Location:\t%d, Size:\t%d, Type:\t%s" % (x[START]+1,x[END]-x[START]+1,x[2]), original_intervals)) 
	log_file.close()

	#print out the genome with snps
	per_line=79
	out_file=open(out_filename,'w')
	print >> out_file, ref_name+"_withINDELS"
	index=0
	while index<len(ref_seq_windels):
		print >> out_file, ref_seq_windels[index:index+per_line]
		index=index+per_line
	out_file.close()	
			
		
