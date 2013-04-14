#!/usr/bin/python


import sys
from copy import deepcopy

#WARNING MAKING ASSUMPTIONS ABOUT EDGES.......
#BAD ASSUMPTIONS

node_to_pos={}
genomic_edges={}
somatic_edges={}
all_edges={}

def to_pos(s):
	s=s.split(':')
	return (to_chr(s[0]),int(s[1]))

def to_chr(s):
	s=s.lower()
	if s[:3]=='chr':
		s=s[3:]
		if s=='x':
			return 23
		if s=='y':
			return 24
	return int(s)

if len(sys.argv)!=4:
	print "%s problem_file paths scale" % sys.argv[0]
	sys.exit(1)



def read_problem_file(filename):
	h=open(filename,'r')
	for line in h:
		line=line.split()
		if line[1]=='NODE':
			#node line
			#c NODE 24 chr1:202709642
			n=int(line[2])
			node_to_pos[n]=to_pos(line[3])
		elif line[1]=='Somatic':
			#Somatic line
			#c Somatic       chr1:202869943  chr1:204056467  2
			f=to_pos(line[2])
			t=to_pos(line[3])
			if f>t:
				print >> sys.stderr, "ERROR" 
				sys.exit(1)
			typ=int(line[4])
			cost=int(line[5])
			if (f,t) not in somatic_edges:
				somatic_edges[(f,t)]=[]
			somatic_edges[(f,t)].append((cost,False,typ))
			somatic_edges[(f,t)].append((cost,False,typ))
		elif line[1]=='Genomic':
			#genomid line
			#c Genomic       chr2:8917257    chr2:8921728
			f=to_pos(line[2])
			t=to_pos(line[3])
			if f>t:
				print >> sys.stderr, "ERROR" 
				sys.exit(1)
			cost=int(line[4])
			if (f,t) not in genomic_edges:
				genomic_edges[(f,t)]=[]
			genomic_edges[(f,t)].append((cost,True,0))
			genomic_edges[(f,t)].append((cost,True,0))
	h.close()
	for (f,t) in somatic_edges:
		if (f,t) not in all_edges:
			all_edges[(f,t)]=[]
		all_edges[(f,t)]+=somatic_edges[(f,t)]
		all_edges[(f,t)].sort()
	for (f,t) in genomic_edges:
		if (f,t) not in all_edges:
			all_edges[(f,t)]=[]
		all_edges[(f,t)]+=genomic_edges[(f,t)]
		all_edges[(f,t)].sort()
	
def annotate_path(p):	
	used={}
	path=[]
	for x in range(1,len(p)):
		f=node_to_pos[p[x-1]]
		t=node_to_pos[p[x]]
		# not the source and sink chromosomes
		if f[0]>0 and t[0]>0:
			#get the canonical
			path.append((f,t))
			c=(f,t)
			if f>t:
				c=(t,f)
			if c not in all_edges:
				print >> sys.stderr, " Dropping edge " ,f,t
				sys.exit(1)
			else:
				if c not in used:
					used[c]={'genomic':0,'somatic':0}
				#find out if we are using somatic or genomic
				genomic=all_edges[c][used[c]['genomic']+used[c]['somatic']][1]	
				#print len(all_edges[c]),used[c],c
				if genomic:
					used[c]['genomic']+=1
				else:
					used[c]['somatic']+=1
	return used,path

def read_paths_file(filename):
	p=[]
	h=open(filename,'r')
	for line in h:
		if line.find('Search ... ')==0:
			p=eval(line.replace('Search ... ',''))
	h.close()
	if len(p)!=0:
		#break up the paths
		sub_paths=[]
		sub_path=[]
		for x in p:
			if x<=6:
				if len(sub_path)>0:
					sub_paths.append(sub_path)
					sub_path=[]
			else:
				sub_path.append(x)
		#annotate each path
		for sub_path in sub_paths:
			print "# " + filename + " subpath_len" , len(sub_path)
			used,path=annotate_path(sub_path)
			l_genomic_edges=[]
			for f,t in genomic_edges:
				u=0
				if (f,t) in used:
					u=used[(f,t)]['genomic']*scale
				l_genomic_edges.append((f,t,0,u))	
			print l_genomic_edges
			l_somatic_edges=[]
			for f,t in somatic_edges:
				u=0
				if (f,t) in used:
					u=used[(f,t)]['somatic']*scale
				l_somatic_edges.append((f,t,somatic_edges[(f,t)][0][2],u))
			print l_somatic_edges
	return p

problem_filename=sys.argv[1]
paths_filename=sys.argv[2]
scale=float(sys.argv[3])

read_problem_file(problem_filename)
read_paths_file(paths_filename)

