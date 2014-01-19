import sys

from random import choice

if len(sys.argv)!=4:
	print "%s flow_file out alpha"
	sys.exit(1)

flow_filename=sys.argv[1]
out_filename=sys.argv[2]
alpha=float(sys.argv[3])

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
		e[fn].pop(tn)
	if len(e[fn])==0:
		e.pop(fn)



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


def print_seq(ls):
	if len(ls)==0:
		return ""
	os=[]
	for s in ls:	
		if len(s)==1:
			os.append(str[0])
			continue
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
		os.append(o)
	return "|".join(os)
		

def is_primary(l):
	l=canonical_loop(l)
	try:
		sa=l.index(1)
		return True
	except:
		pass
	return False	

def loop_to_loops(l):
	l=canonical_loop(l)
	sa=-1
	sb=-1
	#l cannot have repeat element, but can go through source/sink on other side
	try:
		sa=l.index(1)
	except ValueError:
		pass
	try:
		sb=l.index(4)
	except ValueError:
		pass
	if sa>=0 and sb>=0:
		#this is a double loop
		l1=l[:sb]
		l2=canonical_loop(l[sb:])
		return [(tuple(l1),),(tuple(l2),)]
	else:
		return [(tuple(l),)]

all_loops={}


def join_loops(l1,l2):
	nl=set()
	for l in l1:
		nl.add(l)
	for l in l2:
		nl.add(l)
	nl=list(nl)
	nl.sort()
	return tuple(nl)

def loop_to_loop(l1,l2):
	if len(l2)!=1:
		print >> sys.stderr, "FAILEDS AT THIS"
		sys.exit(1)
	l2=l2[0]
	#check if the loop is already in l1
	for loop in l1:
		if loop==l2:
			return ()
	#check for overlaping node
	for loop in l1:
		for x in loop:
			if x in l2:
				#they overlap
				l3=list(l1)+[l2]
				l3.sort()
				return tuple(l3)
	return ()
	
		

#take in contigs
# { '(1,3,4,16,347,17,547,245,724572)':count ...}
#same for loops
def cplexout(lines,loops,candidates,filename):
	#should really remove length 1 and copy count 1 candidates
	new_candidates=[]
	for c in candidates:
		r=[]
		q=[]
		for x in c:
			r.append(x)
			q.append(c[x])
		r=tuple(r)
		q=tuple(q)
		new_candidates.append((r,q))
	candidates=new_candidates
	
	of=open(filename,'w')
	#index the required and candidates
	idxs={}
	line_idxs=[]
	idx_to_line={}
	loop_idxs=[]
	idx_to_loop={}
	candidate_idxs=[]
	idx_to_candidate={}
	for x in lines:
		line_idxs.append(len(idxs))
		idx_to_line[len(idxs)]=x
		idxs[x]=len(idxs)
	for x in loops:
		loop_idxs.append(len(idxs))
		idx_to_loop[len(idxs)]=x
		idxs[x]=len(idxs)
	for x in candidates:
		candidate_idxs.append(len(idxs))
		idx_to_candidate[len(idxs)]=x
		idxs[x]=len(idxs)
	ks=idxs.keys()
	#generate the objective 
	obj=[]
	for x in range(len(idxs)):
		obj.append(str(alpha)+' m'+str(x))
		obj.append('i'+str(x))
	#for x in idxs:
	#	obj.append(str(lmbda)+' m'+str(idxs[x]))
	#for x in idxs:
	#	obj.append('i'+str(idxs[x]))
	print >> of, "Minimize\n obj:"," + ".join(obj)
	
	#print the conditions
	print >> of, "Subject To"
	for x in line_idxs:
		l=idx_to_line[x]
		terms=[]
		terms.append("m"+str(x))
		for y in candidate_idxs:
			c,cp=idx_to_candidate[y]
			if l[0] in c:
				# print M * (how many times c has l in it)
				i=c.index(l[0])
				terms.append(str(cp[i]) + " m"+str(y))
		print >> of,  " l"+str(idxs[l])+":\t" + " + ".join(terms) + " = " + str(lines[l])
	for x in loop_idxs:
		l=idx_to_loop[x]
		terms=[]
		terms.append("m"+str(x))
		for y in candidate_idxs:
			c,cp=idx_to_candidate[y]
			if l[0] in c:
				# print M * (how many times c has l in it)
				i=c.index(l[0])
				terms.append(str(cp[i]) + " m"+str(y))
		print >> of, " o"+str(idxs[l])+":\t" + " + ".join(terms) + " = " + str(loops[l])
	for x in range(len(idxs)):
		if x in idx_to_line:
			print >> of , "\\ LINE " , print_seq(idx_to_line[x])	
		elif x in idx_to_loop:
			print >> of , "\\ LOOP " , print_seq(idx_to_loop[x])	
		elif x in idx_to_candidate:
			print >> of , "\\ CANDIDATE " , print_seq(idx_to_candidate[x][0]), idx_to_candidate[x][1]	
		else:
			print "MAJOR ERROR Xsadf"
			sys.exit(1)
		print >> of,  " lx"+str(x)+":\t" + " m"+ str(x) +" - 1000 i"+str(x)+ " < 0"
		

	#print the bounds
	print >> of, "Bounds"
	for b in range(len(idxs)):
		print >> of, "0 <= m" + str(b) + " <= 1000"

	print >> of, "Binary"
	print >> of, "\t"+ " ".join(map(lambda x : "i"+str(x),range(len(idxs))))
		
	print >> of, "General"
	print >> of, "\t"+ " ".join(map(lambda x : "m"+str(x),range(len(idxs))))

	print >> of , "End"
	of.close()



def common_factor(depth,s):
	for x in range(2,depth):
		good=True
		for loop in s:
			if s[loop]%x!=0:
				good=False
				break
		if good:
			return True
	return False
	

def sloop_candidates(depth,s):
	candidates=sloop_candidates_r(depth,s)
	new_candidates=[]
	for x in candidates:
		if common_factor(depth,x):
			continue
		new_candidates.append(x)
	return candidates

def sloop_candidates_r(depth,s):
	r=[]
	if len(s)==1:
		for x in range(depth):
			r.append({s[0]:(x+1)})
	else:
		qq=sloop_candidates(depth,s[1:])
		for q in qq:
			for x in range(depth):
				d=q.copy()
				d[s[0]]=x+1
				r.append(d)
	return r


def decompose(oe):
	e=edges_copy(oe)
	sz=0
	for fn in e:
		for tn in e[fn]:
			sz+=e[fn][tn]	
	lines={}
	loops={}
	
	#find the neighbours
	stack=[]
	while True:
		if len(stack)<2:
			start_edge=get_start_edge(e)
			if start_edge==False:
				break
			decrease(e,start_edge[0],start_edge[1])
			sz-=1
			stack=[start_edge[0],start_edge[1]]
		fn=stack[-1]
		if len(e[fn])==1:
			tn=e[fn].keys()[0]
			decrease(e,fn,tn)
		else:
			ns=[]
			for tn in e[fn]:
				for x in range(e[fn][tn]):
					ns.append(tn)
			if len(ns)==0:
				break
			#choose a neighbour 
			tn=choice(ns)
			decrease(e,fn,tn)
		sz-=1
		if tn in stack:
			#remove the loop!
			idx=stack.index(tn)
			#loop=canonical_loop(stack[idx:])
			new_loops=loop_to_loops(stack[idx:])
			for loop in new_loops:
				print print_seq(loop),
				d=None
				if is_primary(loop[0]):
					print "LINE"
					d=lines
				else:
					print "LOOP"
					d=loops
				if loop not in d:
					d[loop]=0
					d[(reverse_path(loop[0]),)]=0
				d[loop]+=1
				d[(reverse_path(loop[0]),)]+=1
			stack=stack[:idx]
			#if len(loops)>3 and len(lines)>3:
			#	break
		stack.append(tn)



	#find loop links	
	sloops=set(loops.keys())
	quick_loop_links=set()
	checked_loops=set()
	loop_links=set()
	for loop in sloops:
		for loopy in sloops:
			new_loop=join_loops(loop,loopy)
			if new_loop in checked_loops:
				continue
			checked_loops.add(new_loop)
			if len(loop_to_loop(loop,loopy))>0:
			
				quick_loop_links.add(sum(map(sum,loop))+sum(map(sum,loopy)))
				loop_links.add(new_loop)
	
	#find the line links
	quick_line_loop_links=set()
	line_loop_links=set()
	for line in lines:
		for loop in sloops:
			new_loop=join_loops(line,loop)
			if new_loop in checked_loops:
				continue
			checked_loops.add(new_loop)
			if len(loop_to_loop(line,loop))>0:
				line_loop_links.add(new_loop)
				quick_line_loop_links.add(sum(line[0])+sum(map(sum,loop)))
	
	for x in line_loop_links:
		print "XY",print_seq(x)

	
	#find super loops
	checked_loops=set()
	while True:
		mx=0
		new_loops=[]
		for loop in sloops:
			if loop in checked_loops:
				continue
			checked_loops.add(loop)
			for loopy in sloops:
				if len(loopy)!=1:
					continue
				if loopy[0] in loop:
					continue
				for subloop in loop:
					if loops[loopy]>3*len(loop):
						s=sum(subloop)+sum(loopy[0])
						if s not in quick_loop_links:
							continue
						new_loop=join_loops((subloop,),loopy)
						if new_loop in loop_links:
							#add it
							new_loops.append(join_loops(loop,loopy))
		for new_loop in new_loops:
			mx=max(mx,len(new_loop))
			sloops.add(new_loop)
		print len(sloops),mx
		if len(new_loops)==0:
			break
	
	#add super loops to the candidates
	print "CANDIDATES"
	candidates=[]
	for sloop in sloops:
		print "process candidates..."
		new_candidates=sloop_candidates(3,sloop)
		candidates+=new_candidates
		#check for line attachments
		for line in lines:
			for loop in sloop:
				new_loop=join_loops(line,loop)
				if new_loop in line_loop_links:
					for nc in new_candidates:
						d=nc.copy()
						d[line[0]]=1
						candidates.append(d)
					break
			
					
	#try to print the basic lp solve function
	cplexout(lines,loops,candidates,out_filename)	
	sys.exit(1)	
	




	#find line to loop links
	

	
	#find line to superloops

	#compute candidate contigs

	#print out IP/LP

	#solve

	sys.exit(1)	

	#lets find out which nodes bring flow into secondary loops
	#first lets build dictionary of secondary_loop -> dict
	loops_ins={}
	for loop in loops:
		loops_ins[loop]={}
	for fn in oe:
		for tn in oe[fn]:
			for loop in loops_ins:
				if fn not in loop and tn in loop:
					loops_ins[loop][tn]=oe[fn][tn]

	for loop in loops_ins:
		print print_seq(loop)
		for tn in loops_ins[loop]:
			print "\t"+str(loops_ins[loop][tn])+"\t"+str(tn)
	
	for loop in loops_ins:
		print print_seq(loop)
		c={}
		try:
			for x in range(len(loop)-1):
				y=oe[loop[x]][loop[x+1]]
				if y not in c:
					c[y]=0
				c[y]+=1
			print c
		except:
			loop=reverse_path(loop)
			for x in range(len(loop)-1):
				y=oe[loop[x]][loop[x+1]]
				if y not in c:
					c[y]=0
				c[y]+=1
			print c
			
		
	
				
					
	print sz
	sys.exit(1)
	return loops

def reachable_without(oe,fn,tn):
	e=edges_copy(oe)
	decrease(e,fn,tn)
	visited=set([1,2])
	if tn in (1,2):
		return True
	q=[1,2]
	while len(q)>0:
		n=q.pop() #DFS if pop from back, BFS if pop from front?
		for t in e[n]:
			if tn==t:
				return True
			if t not in visited:
				visited.add(t)
				q.append(t)
	return False
	
		
	

def connected_components_without(oe,fn,tn):
	e=edges_copy(oe)
	decrease(e,fn,tn)
	return unsafe_connected_components(e)

def unsafe_connected_components(e):
	cc=0
	while len(e)>0:
		#grab the first node
		stack=[e.keys()[0]]
		while len(stack)>0:
			fn=stack.pop()
			if fn in e:
				for tn in e[fn]:
					stack.append(tn)
				e.pop(fn)
		cc+=1
	return cc

def connected_components(oe):
	e=edges_copy(oe)
	return unsafe_connected_components(e)

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
