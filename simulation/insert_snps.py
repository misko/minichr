#!/usr/bin/python

import sys
import random

def usage():
	print "%s fasta_reference_filename num_snps output_fasta_file output_log_file" % sys.argv[0]



letters=['A','C','G','T']


if __name__=='__main__':
	if len(sys.argv)!=5:
		usage()
		sys.exit(1)
	
	in_filename=sys.argv[1]
	num_snps=int(sys.argv[2])
	out_filename=sys.argv[3]
	log_filename=sys.argv[4]


	#read in the reference
	in_file=open(in_filename,'rU')
	ref_name=in_file.readline().strip()
	assert(ref_name[0]=='>')
	lines=in_file.readlines()
	ref_seq="".join(map(lambda x : x.strip(), lines))	
	in_file.close()

	#think of some snps
	snp_locations=[]
	snp_base=[]
	assert(num_snps<len(ref_seq))
	r=random.Random()
	j=0
	while j<num_snps:
		candidate_pos=r.randrange(len(ref_seq))
		if not candidate_pos in snp_locations and ref_seq[candidate_pos]==ref_seq[candidate_pos].upper():
			snp_locations.append(candidate_pos)
			j=j+1
			#generate a new base
			x=r.randrange(3)+1
			snp_base.append(letters[(letters.index(ref_seq[candidate_pos].upper())+x)%4].lower())

	#print out the log
	log_file=open(log_filename,'w')
	to_write=[]
	for s in range(num_snps):
		to_write.append([snp_locations[s],"SNP: %d -> %s" % (snp_locations[s]+1,snp_base[s])])
	to_write.sort()
	print >> log_file, "\n".join(map( lambda (x,y) : y, to_write)) 
	log_file.close()

	#print out the genome with snps
	ref_seq_wsnps=ref_seq[:]
	for s in range(num_snps):
		ref_seq_wsnps=ref_seq_wsnps[:snp_locations[s]]+snp_base[s].lower()+ref_seq_wsnps[snp_locations[s]+1:]
	per_line=79
	out_file=open(out_filename,'w')
	print >> out_file, ref_name+"_withSNPS"
	index=0
	while index<len(ref_seq):
		print >> out_file, ref_seq_wsnps[index:index+per_line]
		index=index+per_line
	out_file.close()	
			
		
