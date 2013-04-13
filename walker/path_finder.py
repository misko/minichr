#!/usr/bin/python

import os
import sys
from copy import deepcopy
import random
from multiprocessing import Pool

params={'num_nodes':0,'k':0,'cost':0,'mins':0}
edges={}


#provide extra edge to location from node 3
def print_graph(filename,flow_restriction={},extra=-1):
	sx=0
	arc_lines=[]
	for f in edges:
		for t in edges[f]:
			edges[f][t].sort()
			flowr=len(edges[f][t])
			if f!=1:
				if f in flow_restriction and t in flow_restriction[f]:
					if flow_restriction[f][t]<flowr:
						s=sum(map(lambda x : x[0],edges[f][t][flow_restriction[f][t]:flowr]))
						if f==1:
							print s
						#print >> sys.stderr, "Removing " , flowr-flow_restriction[f][t], " capacity from " , f , " to ", t , " - value " , s
						sx+=s
					flowr=min(flowr,flow_restriction[f][t])
				for x in range(flowr):
					score,low,cap = edges[f][t][x]
					arc_lines.append( "a\t%d\t%d\t%d\t%d\t%d" % (f,t,low,cap,score))
			else:
				if len(edges[f][t])!=1:
					print >> sys.stderr, "Failed to process"
					sys.exit(1)
				score,low,cap = edges[f][t][0]
				if f in flow_restriction and t in flow_restriction[f]:
					#print flow_restriction[f][t]
					low=min(low,flow_restriction[f][t])
					cap=min(cap,flow_restriction[f][t])
				arc_lines.append( "a\t%d\t%d\t%d\t%d\t%d" % (f,t,low,cap,score))
				
	if extra>0:
		arc_lines.append( "a\t%d\t%d\t%d\t%d\t%d" % (3,extra,1,1,0))
	h = open(filename,'w')
	print >> h, "p\tmin\t%d\t%d" % (params['num_nodes'],len(arc_lines))
	for x in range(params['num_nodes']):
		print >> h, "n\t%d\t0" % (x+1)
	for arc_line in arc_lines:
		print >> h, arc_line
	h.close()
	return sx

def get_flow(problem_filename):
	fld={}
	cost=1337
	stream = os.popen("cat " + problem_filename + " | cs2.exe")
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
			edges[f][t].append((score,low,cap))
		elif line[0]=='p':
			#handle program parameters
			line=line.split()
			params['num_nodes']=int(line[2])
			params['num_arcs']=int(line[3])
	return walks

def get_candidates(last_node,flow):
	r=[]
	#want to return edges leaving the last node that have positive flow
	if last_node in flow:
		for to_node in flow[last_node]:
			if flow[last_node][to_node]>0:
				r.append(to_node)
	#print "Finding candidates for " , last_node,r
	return r



def setup():
	#lets get the original cost
	print_graph('out')
	flow,cost=get_flow('out')
	params['cost']=cost

	#remove unused edges
	unused_cost=print_graph('out',flow)
	flow,cost=get_flow('out')
	if cost!=params['cost']:
		cerr << "FATAL"
		sys.exit(1)
	unused_cost=print_graph('out',flow)
	flow,cost=get_flow('out')
	if cost!=params['cost']:
		cerr << "FATAL"
		sys.exit(1)
	params['unused_cost']=unused_cost

	print >> sys.stderr, "Moving along..."
	print "Original cost: ", params['cost']
	return flow



def search(t):
	l,fl,path=t
	#fname='/dev/shm/out_t'+str(os.getpid())
	fname='out'
	#ouc=print_graph(fname+"Y",fl,extra=path[-1])
	#ofl,ocost=get_flow(fname+"Y")
	
	results=[]
	#keep going until we dont hit the sink!
	last_node = path[-1]
	if (last_node==5 or last_node==6) and params['walks']==(path.count(5)+path.count(6)):
		unused_cost=print_graph(fname,fl,extra=last_node)
		local_fl,cost=get_flow(fname)
		print "DONE!",params['walks'],params['mins'],params['cost'],cost,unused_cost,len(path)
		print "Search ... " , path[-10:]
		if cost==0:
			params['mins']+=1
		return results
	candidates=get_candidates(last_node,fl)
	for candidate in candidates:
		flr=deepcopy(fl)
		#some quick sanity checks
		if flr[last_node][candidate]<=0:
			print >> sys.stderr, "ERasdfROR FDSFasfd" 
			sys.exit(1)
		flr[last_node][candidate]-=1
		#lets find out if taking it out gives the same cost
		unused_cost=print_graph(fname,flr,extra=candidate)
		#print path,candidate
		local_fl,cost=get_flow(fname)
		if last_node==3 or (len(path)>2 and path[-2]==3) or (len(path)>3 and path[-3]==3):
			d=unused_cost-params['unused_cost']
			print candidate,"cost:",cost,"cost+uc-param",cost+unused_cost-params['unused_cost'],"uc",unused_cost,params['cost']
			print flr[last_node][candidate]
			print path[-10:]
		if cost<=0:
			d=unused_cost-params['unused_cost']
			#print cost,d,cost+d,params['cost']
			if cost+d==params['cost']:
				#if last_node==1:
				#	print "passed!"
				local_path=deepcopy(path)
				local_path.append(candidate)
				results.append((len(local_path),flr,local_path))
			#if last_node==1:
			#	print "decided!"
			#	sys.exit(1)
		else:
			print >> sys.stderr, "ERROR !!!" , cost
			sys.exit(1)
	return results
	#for result in results:
	#	local_flr=deepcopy(fl)
	#	local_flr[last_node][result]-=1
	#	#print flr,path
	#	local_path=deepcopy(path)
	#	local_path.append(result)
	#	search(local_flr,local_path)
	

params['walks']=read_graph(sys.stdin)

tflow=setup()

#p = Pool(32)

q=[(1,tflow,[3])]
while len(q)>0:
	#q.sort()
	#qtop=[]
	#for x in range(min(32*3,len(q))):
	#	i=int((0.7+0.3*random.random()) * len(q))	
	#	qtop.append(q.pop(i))
	#tresults=p.map(search,qtop)
	#for rt in tresults:
	#	for result in rt:
	#		q.append(result)
	#
	#i=int(random.random() * len(q))
	tpl=q.pop()
	for result in search(tpl):
		q.append(result)
