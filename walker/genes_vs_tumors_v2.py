


cancers={}
patient={}
lookup={}

from multiprocessing import Pool
bins=[72,30,15,6]

cols=4

def get_set(l):
	s=set()
	l=l.replace('-','').split(',')
	for g in l:
		if len(g)>0:
			s.add(g)
	return s

import sys

genes_vs_tumors=[]
tumors_vs_genes=[]
for x in range(cols):
	genes_vs_tumors.append({})
	tumors_vs_genes.append({})

#TCGA-ER-A197    SKCM    gfloopamps      -       -       -       -
#populate per patient per column
for line in sys.stdin:
	line=line.strip().split()
	patient_id=line[0]
	cancer=line[1]
	#if new patient initialize
	if patient_id not in patient:
		patient[patient_id]=[]
		for x in range(cols):
			patient[patient_id].append(set())
	#for each thing update
	for x in range(cols):
		patient[patient_id][x].update(get_set(line[3+x]))
	#update genes_vs_tumors
	for x in range(cols):
		for gene in get_set(line[3+x]):
			if gene not in genes_vs_tumors[x]:
				genes_vs_tumors[x][gene]={}
			if cancer not in genes_vs_tumors[x][gene]:
				genes_vs_tumors[x][gene][cancer]=0
			genes_vs_tumors[x][gene][cancer]+=1
	#update tumors_vs_genes
	for x in range(cols):
		if cancer not in tumors_vs_genes[x]:
			tumors_vs_genes[x][cancer]={}
		for gene in get_set(line[3+x]):
			if gene not in tumors_vs_genes[x][cancer]:
				tumors_vs_genes[x][cancer][gene]=0
			tumors_vs_genes[x][cancer][gene]+=1
	# add to lookup
	if cancer not in cancers:
		cancers[cancer]=set()
	cancers[cancer].add(patient_id)		
	if 'all' not in cancers:
		cancers['all']=set()
	cancers['all'].add(patient_id)		


#trim
trim=3
for x in range(cols):
	rt=set()
	for tumor in tumors_vs_genes[x]:
		v=set()
		for gene in tumors_vs_genes[x][tumor]:
			v.add(tumors_vs_genes[x][tumor][gene])
		if len(v)>2:
			v.remove(max(v))
			v.remove(max(v))
			etrim=max(trim,max(v))
		else:
			etrim=trim
		#etrim=min(etrim,max(len(cancers[tumor])/3,2))
		r=set()
		for gene in tumors_vs_genes[x][tumor]:
			if tumors_vs_genes[x][tumor][gene]<etrim:
				r.add(gene)
		for z in r:
			tumors_vs_genes[x][tumor].pop(z)
		if len(tumors_vs_genes[x][tumor])==0:
			rt.add(tumor)
	for t in rt:
		tumors_vs_genes[x].pop(t)

for x in range(cols):
	rg=set()
	for gene in genes_vs_tumors[x]:
		v=set()
		for tumor in genes_vs_tumors[x][gene]:
			v.add(genes_vs_tumors[x][gene][tumor])
		if len(v)>2:
			v.remove(max(v))
			v.remove(max(v))
			etrim=max(trim,max(v))
		else:
			etrim=trim
		r=set()
		for tumor in genes_vs_tumors[x][gene]:
			if genes_vs_tumors[x][gene][tumor]<etrim:
				r.add(tumor)
		for z in r:
			genes_vs_tumors[x][gene].pop(z)
		if len(genes_vs_tumors[x][gene])==0:
			rg.add(gene)
	for g in rg:
		genes_vs_tumors[x].pop(g)

def simplify(l):
	d={}
	ps={}
	for v,p,x in l:
		if v not in d:
			d[v]=[]
			ps[v]=[]
		d[v].append(x)
		ps[v].append(p)
	k=d.keys()
	k.sort(reverse=True)
	s=[]
	for z in k:
		s.append("("+str(z)+") "+"%0.1e " % ps[z][0]+",".join(d[z]))
	return "\t".join(s)

def subsets(s,sz):
        r=[]
        for e in s:
                if sz==1:
                        r.append(set((e,)))
                else:   
                        ss=subsets(s-set((e,)),sz-1)
                        for sse in ss:
                                r.append(sse.union(set((e,))))
        return r

from random import sample
from math import log 

import math

def nCr(n,r):
    f = math.factorial
    return f(n) / f(r) / f(n-r)

def random_estimate(patients,col,samples,j):
	sets=map(lambda x : sample(patients,j),range(samples))
	#find out the average prob for this class
	sm=0
	for p in patients:
		sm+=len(patient[p][col])/24000.0
	a=sm/len(patients)
	sm=0
	for s in sets:
		#compute the probability
		pro=1.0
		for p in patients:
			z=max(len(patient[p][col])/24000.0,a)
			if p in s:
				pro*=z
			else:
				pro*=(1-z)
		sm+=pro
	return sm/samples


def random_estimate_mp_helper(x):
	patients,col,samples,j,processes=x
	return random_estimate(patients,col,samples,j)
	
from time import sleep

pool = Pool(processes=30) 
def random_estimate_mp(patients,col,samples,j,processes):
	r=sum(pool.map(random_estimate_mp_helper, processes*[(patients,col,samples/processes,j,processes)] ))/processes
	return r

def compute_stats(tumor,col):
	patients=cancers[tumor]
	ps=[]
	for j in range(len(patients)):
		#find out how 
		if j==0:
			ps.append(0)
			continue
		else:
			if j>3 and ps[-1]<=1e-70:
				ps.append(1e-70)
			else:
				ps.append(max(random_estimate_mp(patients,col,1000,j,10)*nCr(len(patients),j),1e-70))
	#sum the tails
	psx=[]
	for j in range(len(patients)):
		psx.append(sum(ps[j:]))
	return psx
	
			

#print tumors_vs_genes
genes_used={}
probs={}
for x in range(cols):
	if x not in probs:
		probs[x]={}
	print "###################"
	print "#TUMORS_VS_GENES_AMP_"+str(bins[x])
	print "###################"
	#first get all genes 
	genes=set()
	for tumor in tumors_vs_genes[x]:
		genes.update(set(tumors_vs_genes[x][tumor].keys()))
	genes=list(genes)
	genes.sort()
	print "Tumor\tGenes"
	for tumor in tumors_vs_genes[x]:
		ps=compute_stats(tumor,x)
		probs[x][tumor]=ps
		print tumor+"("+str(len(cancers[tumor]))+")\t",
		l_genes=[]
		for gene in genes:
			if gene in tumors_vs_genes[x][tumor] and ps[tumors_vs_genes[x][tumor][gene]]<1e-5:
				l_genes.append((tumors_vs_genes[x][tumor][gene],ps[tumors_vs_genes[x][tumor][gene]],gene))
		l_genes.sort(reverse=True)
		if x not in genes_used:
			genes_used[x]=set()
		genes_used[x].update(set(map( lambda x : x[2], l_genes)))
		print simplify(l_genes)
sys.exit(1)
#print genes_vs_tumors
for x in range(cols):
	print "###################"
	print "#GENES_VS_TUMORS_AMP_"+str(bins[x])
	print "###################"
	#first get all tumors
	tumors=set()
	for gene in genes_vs_tumors[x]:
		tumors.update(set(genes_vs_tumors[x][gene].keys()))
	tumors=list(tumors)
	tumors.sort()
	print "Genes\tTumors"
	for gene in genes_vs_tumors[x]:
		if gene in genes_used[x]:
			l_tumor=[]
			for tumor in tumors:
				if tumor in genes_vs_tumors[x][gene] and probs[x][tumor][genes_vs_tumors[x][gene][tumor]]<1e-5:
					l_tumor.append((genes_vs_tumors[x][gene][tumor],probs[x][tumor][genes_vs_tumors[x][gene][tumor]],tumor))
			l_tumor.sort(reverse=True)
			if len(l_tumor)>0:
				print gene+"\t",
				print simplify(l_tumor)
