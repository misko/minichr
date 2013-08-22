#!/usr/bin/python

import os
import sys
from copy import deepcopy
import random
#from multiprocessing import Pool
from subprocess import Popen, PIPE, STDOUT

import time

params={'num_nodes':0,'k':0,'cost':0,'mins':0}
edges={}

pid=os.getpid()

#provide extra edge to location from node 3
def print_graph(flow_restriction={},start_node=None,extra=-1):
	sx=0
	arc_lines=[]
	for f in edges:
		for t in edges[f]:
			score,low,cap=edges[f][t]
			used=0
			if f in flow_restriction and t in flow_restriction[f]:
				used=flow_restriction[f][t]
			if used!=cap:
				low=max(0,low-used)
				arc_lines.append( "a\t%d\t%d\t%d\t%d\t%d" % (f,t,low,cap-used,score))
				sx+=-1*(cap-used)
	if extra>0:
		arc_lines.append( "a\t%d\t%d\t%d\t%d\t%d" % (start_node,extra,1,1,0))
	arc_lines.append("a\t%d\t%d\t%d\t%d\t%d" % (1,1,0,0,0))

	lines=[]
	lines.append("p\tmin\t%d\t%d" % (params['num_nodes'],len(arc_lines)))
	for x in range(params['num_nodes']):
		lines.append("n\t%d\t0" % (x+1))
	for arc_line in arc_lines:
		lines.append(arc_line)
	return sx,"\n".join(lines)

def get_flow(s):
	fld={}
	cost=1337
	#stream = os.popen("cat " + problem_filename + " | /data/misko/2013.04.12/cs2-4.3/cs2.exe")
	
	p = Popen("/filer/misko/cs2.exe",stdin=PIPE,stdout=PIPE)
	stream=p.communicate(input=s)[0].splitlines()
	for line in stream:
		if line.find('cost')>0:
			cost = int(line.split()[4])
		if line[0]=='f':
			line=line.split()
			f=int(line[1])
			t=int(line[2])
			fl=int(line[3])
			if fl>0:
				if f not in fld:
					fld[f]={}
				if t not in fld[f]:
					fld[f][t]=0
				fld[f][t]+=fl
	#print >> sys.stderr, " SOLVED FLOW ", cost
	return fld,cost

def read_graph(x):
	#expect graph
	walks=0
	for line in x:
		if line[0]=='c':
			#comment, essentially skip
			pass
		elif line[0]=='n':
			#handle new node
			pass
		elif line[0]=='a':
			#handle a arc
			line=line.split()
			f=int(line[1])
			t=int(line[2])
			low=int(line[3])
			walks=max(walks,low)
			cap=int(line[4])
			score=int(line[5])
			if f not in edges:
				edges[f]={}
			if t not in edges[f]:
				edges[f][t]=[]
			else:
				"already read this edge"
				sys.exit(1)
			edges[f][t]=(score,low,cap)
		elif line[0]=='p':
			#handle program parameters
			line=line.split()
			params['num_nodes']=int(line[2])
			params['num_arcs']=int(line[3])
	return walks

def get_node(cc, starts, flow):
	nodes=list(rk[cc])
	random.shuffle(nodes)
	for node in nodes:
		if node not in starts:
			for x in edges[node]:
				score,low,cap=edges[node][x]
				used=0
				if node in flow and x in flow[node]:
					used=flow[node][x]
				if used<cap:
					return node
	return None	

def get_candidates(last_node,flow):
	r=[]
	#want to return edges leaving the last node that have positive flow
	for node in edges[last_node]:
		score,low,cap=edges[last_node][node]
		used=0
		if last_node in flow and node in flow[last_node]:
			used=flow[last_node][node]
		if used<cap:
			r.append(node)
	random.shuffle(r)
	return r



max_so_far=[-100000000000]
go_back=0

def search(t):
	fl,path,to_use,starts=t
	if to_use<3000 and to_use>-3000:
		print t
	
	results=[]
	#keep going until we dont hit the sink!
	last_node = path[-1]
	if last_node==starts[-1] and to_use==0:
		unused_cost,g=print_graph(fl,starts[-1],extra=last_node)
		print >> sys.stderr, "DONE!",params['walks'],params['mins'],params['cost'],to_use,unused_cost,len(path)
		if to_use==0:
			r="Search ... " + str(path)
			params['mins']+=1
			return results,True,r
		return results,False,""
	candidates=get_candidates(last_node,fl)
	if len(candidates)==0:
		pass
	elif len(candidates)==1:
		candidate=candidates[0]
		if last_node not in fl:
			fl[last_node]={}
		if candidate not in fl[last_node]:
			fl[last_node][candidate]=0
		fl[last_node][candidate]+=1
		path.append(candidate)
		results.append((fl,path,to_use+edges[last_node][candidate][0],starts))
	else:
		for candidate in candidates:
			flr=deepcopy(fl)
			#some quick sanity checks
			if last_node not in flr:
				flr[last_node]={}
			if candidate not in flr[last_node]:
				flr[last_node][candidate]=0
			flr[last_node][candidate]+=1
				
			#lets find out if taking it out gives the same cost
			max_cost,g=print_graph(flr,starts[-1],extra=candidate)
			#print path,candidate
			local_fl,flow_cost=get_flow(g)
			if max_cost==flow_cost:
				if max_cost>max_so_far[0]+100 or max_cost>-100:
					print max_cost,flow_cost,len(path),to_use
					sys.stdout.flush()
					max_so_far[0]=max_cost
				local_path=deepcopy(path)
				local_path.append(candidate)
				results.append((flr,local_path,to_use+edges[last_node][candidate][0],starts))
			else:
				#print >> sys.stderr, "ERROR !!!" , max_cost,starts
				pass
	return results,False,""
	

params['walks']=read_graph(sys.stdin)
print params
sx,g=print_graph()
flx,cx=get_flow(g)
flx_keys=flx.keys()
random.shuffle(flx_keys)



#try to find the components
k={} #which edge belongs to what component
rk={} #reverse k
cc=0
sz={} #flow per component
for f in flx_keys:
	if f not in k:
		cc+=1
		rk[cc]=set([f])
		sz[cc]=0
		q=[f]
		while len(q)!=0:
			x=q.pop(0)
			new_children=[]
			if x in flx:
				for c in flx[x]:
					if c not in k and flx[x][c]>0:
						k[c]=cc
						rk[cc].add(c)
						new_children.append(c)
						sz[cc]+=flx[x][c]
			q+=new_children

print "There are " ,cc
print sz
#sys.exit(1)	
			


def look(start,end,cc):
	start_node=get_node(cc,[],{})	
	q=[({},[start_node],sz[cc],[start_node])]
	while len(q)>0:
		tpl=q.pop()
		results,d,r = search(tpl)
		if d or time.time()>end:
			return r
		for result in results:
			if result[1][0]==result[1][-1]:
				#print "loop found"
				starts=result[3][:]
				flr=result[0]
				p=result[1]
				to_use=result[2]			
				#now we find a new spot to start from in this cc
				p.append(-10)
				new_start=get_node(cc,starts,flr)
				p.append(new_start)
				starts.append(new_start)
				q.append((result[0],p,to_use,starts))
				print "LOOP"
				q.append(result)
	return ""

s=[]
iz={}
ll={}


def scc(l_edges,v):
	iz[v]=len(iz)
	ll[v]=iz[v]
	s.append(v)
	
	for w in l_edges[v]:
		if w not in iz:
			scc(l_edges,w)
			ll[v]=min(ll[v],ll[w])
		elif w in s:
			ll[v]=min(ll[v],iz[w])
	
	sxx=[]
	if ll[v]==iz[v]:
		while len(s)>0:
			sxx.append(s.pop())
			if sxx[-1]==v:
				break
	return sxx

def tarjan_scc(l_edges):
	nodes=set()
	for f in l_edges:
		nodes.add(f)
		for t in l_edges[f]:
			nodes.add(t)
	sys.setrecursionlimit(2*len(nodes))
	sccs=[]
	for v in nodes:
		if v not in iz:
			sccs.append(scc(l_edges,v))
	return sccs


def cloop(l):
	return tuple(l[l.index(max(l)):]+l[:l.index(max(l))])

def dfs(flr,found):
	nodes=set()
	for f in edges:
		nodes.add(f)
		for t in edges[f]:
			nodes.add(t)
	rnodes=list(nodes)
	random.shuffle(rnodes)
	pi={}
	seen=set()
	for v in rnodes:
		q=[([],v)]
		while len(q)>0:
			parents,node=q.pop()
			for t in edges[node]:
				used=0
				if node in flr and t in flr[node]:
					used=flr[node][t]
				if edges[node][t][2]>used:
					if t in parents:
						#cycle?
						p=cloop(parents[parents.index(t):]+[node])
						if p not in found:
							return p
						else:
							#ignore edge
							pass
					elif t not in seen:
						seen.add(t)
						q.append((parents+[node],t))
	return None

def remove_cycles():
	flr={}
	l_edges=deepcopy(edges)
	taboo=set()
	go_on=True
	
	sxxs=tarjan_scc(l_edges)
	for sxx in sxxs:
		#try to remove largest multiple of this from graph to keep flow
		m=1
		#check if can remove 1, should be able to , sanity check
		for x in range(1,len(sxx)):
			f=sxx[x-1]
			t=sxx[x]
			print f,t
			if t not in edges or f not in edges[t]:
				print "ERROR",edges[t],edges[f]
				sys.exit(1)
	print sxxs


found={}
flr={}
for f in edges:
	flr[f]={}
	for t in edges[f]:
		flr[f][t]=0
k=dfs(flr,found)
while k!=None:
	#try to make new edge list
	tk=tuple(k)
	if tk not in found:
		found[tk]=0
	else:
		print "EERROROROROROOROROROROR"
		sys.exit(1)
	cap_ok=True
	m=0
	while cap_ok:
		flr_p=deepcopy(flr)
		for x in range(len(k)):
			f=k[x]
			t=k[(x+1)%len(k)]
			#make sure that we can subtract this based on cap
			used=flr_p[f][t]
			if edges[f][t][2]>used:
				flr_p[f][t]+=1
			else:
				cap_ok=False
				break
		if cap_ok:
			#lets try to flow solve
			#print "trying to solve"
			max_cost,g=print_graph(flr_p)
			local_fl,flow_cost=get_flow(g)
			if flow_cost==max_cost:
				m+=1
				#print "upping cycle count",m
				flr.update(flr_p)
			else:
				#print "cannot up cycles any more",m
				cap_ok=False
		else:
			#print "no good cap"
			max_cost,g=print_graph(flr)
			local_fl,flow_cost=get_flow(g)
			#print max_cost,flow_cost
		found[tk]=m
	k=dfs(flr,found)

for k in found:
	print len(k),found[k]
sys.exit(1)				

origin_time = time.time()
max_time=60*60*60

d=60*60*12
while True:
	b=[]
	for cc in sz:
		print "working on cc ",cc
		print >> sys.stderr, "working on cc ",cc
		start = time.time()
		r=look(start,start+d,cc)
		#if time.time()<start+d:
		if len(r)>0:
			#found a solution
			print "FOUND A SOLUTION"
			print >> sys.stderr, "FOUND A SOLUTION"
			b.append(r)
		else:
			#didnt find one
			#d=int(d*1.2)
			print "NO SOLUTION FOUND"
			print >> sys.stderr, "NO SOLUTION FOUND"
		print >> sys.stderr, "Restarting with ", d
	if len(b)>0:
		print "\n".join(b)
		sys.stdout.flush()
