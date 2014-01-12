import sys

from random import choice

if len(sys.argv)!=2:
	print "%s flow_file"
	sys.exit(1)

flow_filename=sys.argv[1]


flow_file=open(flow_filename)

edges={}
ins={}
outs={}

nodes=set()
funky_nodes=set()

for line in flow_file:
	line=line.strip().split()
	if line[0]=='f':
		fl=int(line[3])
		if fl==0:
			continue
		fn=int(line[1])
		tn=int(line[2])
		if fn not in edges:
			edges[fn]={}
		if tn not in edges[fn]:
			edges[fn][tn]=0
		edges[fn][tn]+=fl
		if tn not in ins:
			ins[tn]=0
		ins[tn]+=fl
		if fn not in outs:
			outs[fn]=0
		outs[fn]+=fl
		nodes.add(fn)
		nodes.add(tn)
		if abs(fn-tn)!=2:
			funky_nodes.add(fn)

for x in ins:
	if ins[x]!=outs[x]:
		print >> sys.stderr, "FAIL",x
		sys.exit(1)

for x in outs:
	if ins[x]!=outs[x]:
		print >> sys.stderr, "FAIL",x
		sys.exit(1)
		

def get_start_edge(e):
	es=[]
	for fn in e:
		for tn in e[fn]:
			if e[fn][tn]<=0:
				print >> sys.stderr, "FAILDE!"
				sys.exit(1)
			for x in range(e[fn][tn]):
				es.append((fn,tn))
	if len(es)==0:
		return False
	return choice(es)			

def decrease(e,fn,tn):
	e[fn][tn]-=1
	if e[fn][tn]==0:
		e[fn].pop(tn,None)
	if len(e[fn])==0:
		e.pop(fn,None)



def reverse_path(p):
	p=list(p)
	for x in range(len(p)):
		if p[x]%2==0:
			p[x]-=1
		else:
			p[x]+=1
	p.reverse()
	return tuple(p)
	

def canonical_loop(l):
	mn=min(l)
	if mn%2==0:
		l=reverse_path(l)
	l=l[l.index(min(l)):]+l[:l.index(min(l))] #bring min to front
	return tuple(l)


#given loop and subseq s, give a histogram of things to come next
def ocount(l,s):
	c={}
	sz=len(l)
	szz=len(s)
	for i in range(sz):
		missed=False
		for x in range(szz):
			if l[(i+x)%sz]!=s[x]:
				missed=True
				break
		if not missed:
			z=l[(i+szz)%sz]
			if z!=None:
				if z not in c:
					c[z]=0
				c[z]+=1	
	return c

def is_subsequence(s,m):
	sz=len(m)
	szz=len(s)
	for i in range(sz-szz+1):
		missed=False
		for x in range(szz):
			if m[i+x]!=s[x]:
				missed=True
				break
		if not missed:
			return True
	return False


def print_seq(s):
	if len(s)==0:
		return ""
	if len(s)==1:
		return str(s[0])
	o=""
	d=None
	l=0
	for x in range(len(s)-1):
		if s[x]-s[x+1]!=d:
			if len(o)>0:
				if l>0:
					o+="-"+str(s[x])
				else:
					o+=","+str(s[x])
				#o+=","+str(s[x+1])
			else:
				o=str(s[x])
			d=s[x]-s[x+1]
			l=0
		else:
			l+=1
	o+="-"+str(s[-1])
	return o
			

all_loops={}


def decompose(e):	
	loops={}
	#find the neighbours
	stack=[]
	while True:
		if len(stack)<2:
			start_edge=get_start_edge(e)
			if start_edge==False:
				break
			decrease(e,start_edge[0],start_edge[1])
			stack=[start_edge[0],start_edge[1]]
		fn=stack[-1]
		ns=[]
		for tn in e[fn]:
			for x in range(e[fn][tn]):
				ns.append(tn)
		if len(ns)==0:
			break
		#choose a neighbour 
		tn=choice(ns)
		decrease(e,fn,tn)
		if tn in stack:
			#remove the loop!
			idx=stack.index(tn)
			loop=canonical_loop(stack[idx:])
			if loop not in all_loops:
				loops[loop]=0
				loops[reverse_path(loop)]=0
			loops[loop]+=1
			loops[reverse_path(loop)]+=1
			stack=stack[:idx]
		stack.append(tn)
	return loops

def edges_copy(e):
	ne={}
	for x in e:
		ne[x]=e[x].copy()
	return ne


all_loops=decompose(edges_copy(edges))


#get spectrum

d={}
c={}
k=20
for loop in all_loops:
	for x in range(len(loop)+k-1):
		z=[]
		for y in range(k):
			z.append(loop[(x+y)%len(loop)])
		n=z[0]
		nr=n
		if n%2==0:
			nr-=1
		else:
			n+=1
		if n in funky_nodes or nr in funky_nodes:
			z=tuple(z)
			if z not in d:
				d[z]=0
			d[z]+=1
			if z[0] not in c:
				c[z[0]]={}
			if z[-1] not in c[z[0]]:
				c[z[0]][z[-1]]=0
			c[z[0]][z[-1]]+=1

ks=[]
for k in d:
	ks.append((d[k],k))

ks.sort()
print ks


sys.exit(1)
		


#do some funky stuff

finals=set()

cutoff=0.3
for node in funky_nodes:
	
	go_on=True

	#check if already found somewhere
	for f in finals:
		if node in f:
			go_on=False
			break
	if go_on==False:
		continue

			
	dist={(node,):1}
	while go_on:
		go_on=False # need to find one with higher then cutoff to go on
		new_dist={}
		sum=0
		for s in dist:
			for loop in all_loops:	
				oc=ocount(loop,s)
				for x in oc:
					k=s+(x,)
					if k not in new_dist:
						new_dist[k]=0.0
					new_dist[k]+=oc[x]
					sum+=oc[x]*all_loops[loop]
		for s in new_dist:
			new_dist[s]/=sum
			if new_dist[s]>=cutoff:
				go_on=True
				dist=new_dist
	for x in dist:
		if len(x)>2:
			#check if this supers any other finals
			to_remove=set()
			for f in finals:
				if is_subsequence(f,x):
					to_remove.add(f)
			for f in to_remove:
				finals.remove(f)
			finals.add(x)
			print print_seq(x)


for f in finals:
	print print_seq(f)	
