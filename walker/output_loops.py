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
				print >> sys.stderr, " Dropping edge " ,f, "(",p[x-1],")",t,"(",p[x],")"
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


def canonical_loop(l):
	if len(l)==0:
		return l
	mn=min(l)
	ix=l.index(mn)
	lr=tuple(l[ix:]+l[:ix])
	if mn%2==0 or (len(lr)>2 and lr[-1]<lr[1]):
		lr=map(lambda x : x-1 if x%2==0 else x+1, lr)
		lr.reverse()
		return canonical_loop(tuple(lr))
	return lr

def loops_in_path(p):
	loops={}
	p_so_far=[]
	for x in p:
		#not one of the real nodes...
		if x<=6:
			p_so_far.clear()
			#p_so_far=[]
		try:
			ix=p_so_far.index(x)
			#found a copy
			loop=p_so_far[ix+1:]+[x]
			p_so_far=p_so_far[:ix+1]
			#rotate list to min start
			cloop=canonical_loop(loop)
			#cloop=tuple(loop)
			#print cloop,p
			if cloop not in loops:
				loops[cloop]=0
			loops[cloop]+=1
			#loops.add(canonical_loop(loop))
		except ValueError:
			p_so_far.append(x)
	return loops
			

def read_paths_file(filename):
	alld={}
	samples=0
	h=open(filename,'r')
	for line in h:
		if line.find('Search ... ')==0:
			try:
				p=eval(line.replace('Search ... ',''))
				p_loops=loops_in_path(p)
				#print "XX",len(p_loops)
				samples+=1
				for x in p_loops:
					if x not in alld:
						alld[x]=[]
					alld[x].append(p_loops[x])
				#print "\n".join(map(str,x))
			except Exception, e:
				print >> sys.stderr,  "XX:", e
	h.close()
	#for loop in alld:
	#	print loop
	#sys.exit(1)
	if len(alld)!=0:
		for loop in alld:
			if sum(alld[loop])>3:
				used,path=annotate_path(loop)
				l_genomic_edges=[]
				cq=0
				length=0
				for f,t in genomic_edges:
					u=0
					if (f,t) in used:
						cq+=f[1]
						cq+=t[1]
						u=used[(f,t)]['genomic']*scale*(sum(alld[loop])/len(alld[loop]))#*alld[loop]/float(samples)
						l_genomic_edges.append((f,t,0,u))
						length+=abs(f[1]-t[1])
				l_somatic_edges=[]
				for f,t in somatic_edges:
					u=0
					if (f,t) in used:
						u=used[(f,t)]['somatic']*scale*(sum(alld[loop])/len(alld[loop]))
						l_somatic_edges.append((f,t,somatic_edges[(f,t)][0][2],u))
				print "# " + filename + " subpath_len" , len(loop) ,  alld[loop] , "/", samples, len(alld)
				print "#Samples: %d GlobalFreq: %0.2f Expected: %0.2f Expected if present: %0.2f L: %d" % ( samples , float(len(alld[loop]))/float(samples) , float(sum(alld[loop]))/samples , float(sum(alld[loop]))/len(alld[loop]), length), loop
				#print str(map(lambda x : "(%s, %s, %d, %0.2f)" % (str(x[0]),str(x[1]),x[2],x[3]) , l_genomic_edges)).replace('\'','')
				#print str(map(lambda x : "(%s, %s, %d, %0.2f)" % (str(x[0]),str(x[1]),x[2],x[3]) , l_somatic_edges)).replace('\'','')
	return None

problem_filename=sys.argv[1]
paths_filename=sys.argv[2]
scale=float(sys.argv[3])

read_problem_file(problem_filename)
read_paths_file(paths_filename)

