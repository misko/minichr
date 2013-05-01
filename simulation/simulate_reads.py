#!/usr/bin/python
import sys, random

qual_range=10
qual_offset=ord('I')
qual_std=3
qual_err_std=2

def get_quality():
	c=-1
	while c<1 or c>qual_range:
		c=random.gauss(0,qual_std)
		c=int(c)
	c=chr(qual_offset+qual_range-c)
	return c

def get_error_quality():
	c=-1
	while c<1 or c>qual_range:
		c=random.gauss(0,qual_err_std)
		c=int(c)
	c=chr(qual_offset+c-1)
	return c

# Set our error rates as seen experimentally.
def get_error_rates():
	error_rates = []

	error_rates.append(0.02828)
	error_rates.append(0.02828)
	error_rates.append(0.03275)
	error_rates.append(0.03275)
	error_rates.append(0.05993)
	error_rates.append(0.05993)
	error_rates.append(0.04159)
	error_rates.append(0.04159)
	error_rates.append(0.05107)
	error_rates.append(0.05107)
	error_rates.append(0.03598)
	error_rates.append(0.03598)
	error_rates.append(0.03529)
	error_rates.append(0.03529)
	error_rates.append(0.05057)
	error_rates.append(0.05057)
	error_rates.append(0.04657)
	error_rates.append(0.04657)
	error_rates.append(0.04828)
	error_rates.append(0.04828)
	error_rates.append(0.03723)
	error_rates.append(0.03723)
	error_rates.append(0.03903)
	error_rates.append(0.03903)
	error_rates.append(0.05967)
	error_rates.append(0.05967)
	error_rates.append(0.05645)
	error_rates.append(0.05645)
	error_rates.append(0.05509)
	error_rates.append(0.05509)
	error_rates.append(0.04153)
	error_rates.append(0.04153)
	error_rates.append(0.04106)
	error_rates.append(0.04106)
	error_rates.append(0.06811)
	error_rates.append(0.06811)
	error_rates.append(0.06778)
	error_rates.append(0.06778)
	error_rates.append(0.05457)
	error_rates.append(0.05457)
	error_rates.append(0.04691)
	error_rates.append(0.04691)
	error_rates.append(0.04875)
	error_rates.append(0.04875)
	error_rates.append(0.07017)
	error_rates.append(0.07017)
	error_rates.append(0.07718)
	error_rates.append(0.07718)
	error_rates.append(0.05374)
	error_rates.append(0.05374)
	error_rates.append(0.04939)
	error_rates.append(0.04939)
	error_rates.append(0.04958)
	error_rates.append(0.04958)
	error_rates.append(0.06946)
	error_rates.append(0.06946)
	error_rates.append(0.07782)
	error_rates.append(0.07782)

# 29-34 are just 23-28, since SHRiMP does local align and we haven't good stats
#	error_rates.append(error_rates[23])
#	error_rates.append(error_rates[24])
#	error_rates.append(error_rates[25])
#	error_rates.append(error_rates[26])
#	error_rates.append(error_rates[27])
#	error_rates.append(error_rates[28])
#do whatever
	for x in range(read_len-len(error_rates)):
		error_rates.append(error_rates[-1])

	scale_by=error_rate*read_len/(sum(error_rates[:read_len]))
	for x in range(read_len):
		error_rates[x]=error_rates[x]*scale_by
	
	return (error_rates)

#def get_indel_len():
#	runif = random.random()
#
#	if runif < 0.43:
#		return 1
#	elif runif < 0.43 + 0.17:
#		return 2
#	elif runif < 0.43 + 0.17 + 0.09:
#		return 3 	
#	elif runif < 0.43 + 0.17 + 0.09 + 0.07:
#		return 4
#	elif runif < 0.43 + 0.17 + 0.09 + 0.07 + 0.06:
#		return 5
#	elif runif < 0.43 + 0.17 + 0.09 + 0.07 + 0.06 + 0.05:
#		return 6
#	elif runif < 0.43 + 0.17 + 0.09 + 0.07 + 0.06 + 0.05 + 0.045:
#		return 7
#	elif runif < 0.43 + 0.17 + 0.09 + 0.07 + 0.06 + 0.05 + 0.045 + 0.035:
#		return 8
#	elif runif < 0.43 + 0.17 + 0.09 + 0.07 + 0.06 + 0.05 + 0.045 + 0.035 + 0.03:
#		return 9
#	else:
#		return 10

def random_base():
	return ('A', 'C', 'G', 'T')[random.randint(0,3)]

def random_other_base(base):
	set = ()
	if base == 'A':
		set = ('C', 'G', 'T')
	elif base == 'C':
		set = ('A', 'G', 'T')
	elif base == 'G':
		set = ('A', 'C', 'T')
	elif base == 'T': 
		set = ('A', 'C', 'G')
	elif base == 'N':
		set = ('N', 'N', 'N')   # this will be thrown out anyway...
	else:
		print >> sys.stderr, "ERROR: random_other_base: (%c)" % (base)
		sys.exit(1)

	return set[random.randint(0,2)] 

def random_other_colour(base):
	set = ()
	if base == '0':
		set = ('1', '2', '3')
	elif base == '1':
		set = ('0', '2', '3')
	elif base == '2':
		set = ('0', '1', '3')
	elif base == '3': 
		set = ('0', '1', '2')
	else:
		print >> sys.stderr, "ERROR: random_other_base: (%c)" % (base)
		sys.exit(1)

	return set[random.randint(0,2)] 


rev={'A':'T','a':'t','C':'G','c':'g','G':'C','g':'c','T':'A','t':'a','N':'N','n':'n'}

def reverse(s):
        ret=[]
        for x in range(len(s)):
		if rev.has_key(s[-x-1]):
                	ret.append(rev[s[-x-1]])
		else:
			ret.append(s[-x-1])
        return "".join(ret)


# create a permuted read 
def permute_read(genome, offset, wantlen):
	#read = genome[offset : offset + wantlen * 3]
	#if len(read) < wantlen:
	#	return (None,None,None,None)
	read = genome[offset : offset + wantlen]


	#First lets do the indels, then put in the SNPs
	dels = []
	ins = []
	if (MAX_indels>0):
		wanted_indels=random.randrange(MAX_indels)
		for x in range(wanted_indels):
			position=random.randrange(read_len-MAX_indel_size/2)
			size=random.randrange(MAX_indel_size)+1
			if random.random()>0.5:
				#do a insert
				#generate the insert
				insert=[]
				for y in range(size):
					insert.append(random_base())
				#insert it
				read=read[:position]+"".join(insert)+read[position:]
				ins.append(size)
			else:
				#do a deletion
				dels.append(size)
				read=read[:position]+read[position+size:]
			if len(read)<read_len:
				return (None,None,None,None,None)
				
	#time to put in the SNPs
	snps=0
	if (MAX_SNPs>0):
		wanted_snps=random.randrange(MAX_SNPs)
		avaliable_positions=range(read_len)
		random.shuffle(avaliable_positions)
		snp_positions=avaliable_positions[:wanted_snps]
		for position in snp_positions:
			read=read[:position]+random_other_base(read[position])+read[position+1:]
		snps = len(snp_positions)

	new_read=read[0:wantlen]
	if (len(new_read)!=wantlen):
		return (None,None,None,None,None)

	#make quality values
	q=""
	for x in range(wantlen):
		q+=get_quality()
	return (new_read,snps,dels,ins,q)	


def colourise_read(read):
	colourmat = (	(0, 1, 2, 3), (1, 0, 3, 2), (2, 3, 0, 1), (3, 2, 1, 0) )
	lettermap = { "A" : 0, "a" : 0, "C" : 1, "c" : 1, "G" : 2, "g" : 2, "T" : 3, "t" : 3 }
	last = 'T'

	csread = []

	csread.append('T')
	for c in read:
		csread.append(str(colourmat[lettermap[last]][lettermap[c]]))
		last = c

	return "".join(csread)

def errorise_read(ls_read, error_rates, q):
	# add cs errors to a colourised read
	errs = 0
	i = 0
	new_qual = []
	new_ls_read = []
	cs_read = colourise_read(ls_read)
	new_cs_read = [cs_read[0]]
	while i < len(ls_read):
		runif = random.random()
		letter = ls_read[i]
		number = cs_read[i+1]
		quality = q[i]
		#if runif < error_rates[i - 1]:
		#	letter = random_other_base(letter).lower()
		#	number = random_other_colour(number).lower()
		#	quality = get_error_quality()
		#	errs += 1
		new_qual.append(quality)
		new_ls_read.append(letter)
		new_cs_read.append(number)
		i += 1

	return ("".join(new_ls_read), "".join(new_cs_read), errs,"".join(new_qual))
				

c={'0':'A','1':'C','2':'G','3':'T','.':'N','T':''}
def bwa_format(s):
	r=""
	for x in s[2:]:
		r=r+c[x]	
	return r


def generate_reads(genome, num_reads, read_len):
	#get the error rates and tell the user about them
	error_rates = get_error_rates()

	#tell the user what is going on
	s="contig_file: %s\nnum_reads: %s\nread_len: %s\nMAX_indels: %d\tMAX_indel_size: %d\tMAX_SNPs: %d\n" % (contig_file,num_reads,read_len,MAX_indels,MAX_indel_size,MAX_SNPs)
	s=s+"error_rate: %.3f\tinsert_size_mean: %d\tinsert_size_std: %d\n" % (error_rate,insert_size_mean,insert_size_std)
	s+='Position specific error rates:\n'
	s+="["
	for x in error_rates:
		s=s+'%.4f' % x +','
	s=s[:-1]+']\n'

	#Write a copy to the log
	h=open('generate_log','w')
	h.write(s)
	h.close()
	print >> sys.stderr, s,

	genome_len = len(genome)

	#open some pipes
	#btie_out_left=open('btie_sim_left.fq','w')
	#btie_out_right=open('btie_sim_right.fq','w')
	#bwa_out_left=open('bwa_sim_left.fq','w')
	#bwa_out_right=open('bwa_sim_right.fq','w')
	#bfast_out=open('bfast_sim.fq','w')
	#bfast_out_single=open('bfast_sim_single.fq','w')
	#shrimp_out=open('shrimp_sim.fa','w')
	#shrimp_qual=open('shrimp_qual.fa','w')
	#shrimp_out_single=open('shrimp_out_single.fa','w')
	#cs_shrimp_out_fastq=open('cs_reads.fq','w')
	ls_f_shrimp_out_fastq=open('ls_f_reads.fq','w')
	ls_r_shrimp_out_fastq=open('ls_r_reads.fq','w')
	#cs_shrimp_out_fasta=open('cs_reads.fa','w')
	#ls_shrimp_out_fasta=open('ls_reads.fa','w')
	#cs_shrimp_out_qual=open('cs_reads.qual','w')
	#ls_shrimp_out_qual=open('ls_reads.qual','w')
	#reference_out=open('reference_sim.info','w')

	#while there are not enough reads, keep generating
	read_id=0	
	while read_id<num_reads:
		#Generate a mate pair
		#Get the coordinates of the ends
		offset_left = random.randint(0, genome_len - read_len*3);
		insert_size = int(random.normalvariate(insert_size_mean,insert_size_std))
		offset_right = offset_left+insert_size;
		if (offset_right>=(genome_len-read_len*3)):
			continue
		#offset_right = random.randint(0, genome_len - read_len);
		#offset_left = random.randint(0, genome_len - read_len);
	
		#get the left and right original reads
		(read_left, snps_left, dels_left, ins_left, q_left) = permute_read(genome, offset_left, read_len)
		if read_left == None or read_left.find('N') != -1:
			continue
		(read_right, snps_right, dels_right, ins_right, q_right) = permute_read(genome, offset_right, read_len)
		if read_right == None or read_right.find('N') != -1:
			continue
		read_right=reverse(read_right)
		q_right=reverse(q_right)

		#insert errors
		(ls_read_left,cs_read_left,errs_left,q_left) = errorise_read(read_left,error_rates, q_left) 
		(ls_read_right,cs_read_right,errs_right,q_right) = errorise_read(read_right,error_rates, q_right)

		#write out shrimp, paired reads, in fasta
		s="> READ_"+str(read_id)+"_pos"+str(offset_left+1)+":1"+'\n'+ls_read_left+'\n'+"> READ_"+str(read_id)+"_pos"+str(offset_right+1)+":2"+'\n'+ls_read_right+'\n'
		#ls_shrimp_out_fasta.write(s)
		s="> READ_"+str(read_id)+":1"+'\n'+q_left+'\n'+"> READ_"+str(read_id)+":2"+'\n'+q_right+'\n'
		#ls_shrimp_out_qual.write(s)
		s="> READ_"+str(read_id)+"_pos"+str(offset_left+1)+":1"+'\n'+cs_read_left+'\n'+"> READ_"+str(read_id)+"_pos"+str(offset_right+1)+":2"+'\n'+cs_read_right+'\n'
		#cs_shrimp_out_fasta.write(s)
		s="> READ_"+str(read_id)+":1"+'\n'+q_left+'\n'+"> READ_"+str(read_id)+":2"+'\n'+q_right+'\n'
		#cs_shrimp_out_qual.write(s)

		#write out shrimp, paired reads, in fastq
		#s="@READ_"+str(read_id)+"_pos"+str(offset_left+1)+":1"+'\n'+ls_read_left+'\n+\n'+q_left+"\n"
		s="@READ_"+str(read_id)+'\n'+ls_read_left+'\n+\n'+q_left+"\n"
		ls_f_shrimp_out_fastq.write(s)
		#s="@READ_"+str(read_id)+"_pos"+str(offset_right+1)+":2"+'\n'+ls_read_right+'\n+\n'+q_right+'\n'
		s="@READ_"+str(read_id)+'\n'+ls_read_right+'\n+\n'+q_right+'\n'
		ls_r_shrimp_out_fastq.write(s)
		s="@READ_"+str(read_id)+"_pos"+str(offset_left+1)+":1"+'\n'+cs_read_left+'\n+\n'+q_left+"\n@READ_"+str(read_id)+"_pos"+str(offset_right+1)+":2"+'\n'+cs_read_right+'\n+\n'+q_right+'\n'
		#cs_shrimp_out_fastq.write(s)

		#make the tags
		tag = "READ_%d:1 offset=%d snps=%d errs=%d ins=%s dels=%s | " % (read_id,offset_left + 1, snps_left, errs_left, ins_left, dels_left)	
		#make the right read
		tag+= "READ_%d:2 offset=%d snps=%d errs=%d ins=%s dels=%s" % (read_id,offset_right + 1, snps_right, errs_right, ins_right, dels_right)
		s=tag+'\n'
		#write out reference
		#reference_out.write(s)

		read_id=read_id+1
		if num_reads/1000>0 and ( (read_id % (num_reads/1000)==0 or read_id==num_reads)):
			sys.stderr.write("\r%7.2f%c,\tGenerating read %7.d" % (read_id*100/float(num_reads),'%',read_id))
			sys.stderr.flush()
	sys.stderr.write("\n")



#Default settings
num_reads=1000
read_len=100
MAX_SNPs=5
MAX_indels=2
MAX_indel_size=10
error_rate=0.04
insert_size_mean=500
insert_size_std=70


nargs=len(sys.argv)
if nargs < 2 :
	print >> sys.stderr, "usage: %s [contig_file] [num_reads=1000] [read_len=100] [max_indels=2] [max_indel_size=10] [max_snps=5] [error_rate=0.04] [insert_size_mean=500] [insert_size_std=70]" % (sys.argv[0])
	sys.exit(1)

#parse in command line args
contig_file = sys.argv[1]
#optional arguments, i.e. have default values
if nargs>2:
	num_reads = int(sys.argv[2])
if nargs>3:
	read_len = int(sys.argv[3])
if nargs>4:
	MAX_indels = int(sys.argv[4])
if nargs>5:
	MAX_indel_size = int(sys.argv[5])
if nargs>6:
	MAX_SNPs = int(sys.argv[6])
if nargs>7:
	error_rate = float(sys.argv[7])
if nargs>8:
	insert_size_mean = int(sys.argv[8])
if nargs>9:
	insert_size_std = int(sys.argv[9])

print >> sys.stderr, "Loading genome..."

#Read in the genome
fd = open(contig_file, "r")
lines=fd.readlines()
assert(lines[0][0]=='>')
genome = "".join(lines[1:]).replace('\n','').upper()
fd.close()

#inform user and start generating reads
print >> sys.stderr, "Genome loaded (%d bases)." % len(genome)
print >> sys.stderr, "Generating %d reads..." % (num_reads)
generate_reads(genome, num_reads, read_len)
