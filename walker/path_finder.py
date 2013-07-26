#!/usr/bin/python

import os
import sys
from copy import deepcopy
import random
from multiprocessing import Pool

import time

start_node=0
params={'num_nodes':0,'k':0,'cost':0,'mins':0}
edges={}

pid=os.getpid()

#provide extra edge to location from node 3
def print_graph(filename,flow_restriction={},extra=-1):
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
	#stream = os.popen("cat " + problem_filename + " | /data/misko/2013.04.12/cs2-4.3/cs2.exe")
	stream = os.popen("cat " + problem_filename + " | /filer/misko/cs2.exe")
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
	stream.close()
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



def setup():
	#lets get the original cost
	fname='out_%d' % pid
	print_graph(fname)
	flow,cost=get_flow(fname)
	params['cost']=cost
	
	#remove unused edges
	unused_cost=print_graph(fname,flow)
	flow,cost=get_flow(fname)
	if cost!=params['cost']:
		cerr << "FATAL"
		sys.exit(1)
	sys.exit(1)
	unused_cost=print_graph(fname,flow)
	flow,cost=get_flow(fname)
	if cost!=params['cost']:
		cerr << "FATAL"
		sys.exit(1)
	params['unused_cost']=unused_cost

	print >> sys.stderr, "Moving along..."
	print "Original cost: ", params['cost']
	return flow

max_so_far=[-100000000000]

def search(t):
	l,fl,path,cost=t
	#fname='/dev/shm/out_t'+str(os.getpid())
	fname='/dev/shm/out_%d' % pid
	#ouc=print_graph(fname+"Y",fl,extra=path[-1])
	#ofl,ocost=get_flow(fname+"Y")
	
	results=[]
	#keep going until we dont hit the sink!
	last_node = path[-1]
	#if (last_node==5 or last_node==6) and params['walks']==(path.count(5)+path.count(6)):
	#local_fl,cost=get_flow(fname)
	if last_node==start_node and cost==0:
		unused_cost=print_graph(fname,fl,extra=last_node)
		print >> sys.stderr, "DONE!",params['walks'],params['mins'],params['cost'],cost,unused_cost,len(path)
		if cost==0:
			print "Search ... " , path
			params['mins']+=1
			return results,True
		return results,False
	candidates=get_candidates(last_node,fl)
	for candidate in candidates:
		flr=deepcopy(fl)
		#some quick sanity checks
		if last_node not in flr:
			flr[last_node]={}
		if candidate not in flr[last_node]:
			flr[last_node][candidate]=0
		flr[last_node][candidate]+=1
		#lets find out if taking it out gives the same cost
		unused_cost=print_graph(fname,flr,extra=candidate)
		#print path,candidate
		local_fl,cost=get_flow(fname)
		#if last_node==3 or (len(path)>2 and path[-2]==3) or (len(path)>3 and path[-3]==3):
		#	d=unused_cost-params['unused_cost']
		#	#print candidate,"cost:",cost,"cost+uc-param",cost+unused_cost-params['unused_cost'],"uc",unused_cost,params['cost']
		#	#print flr[last_node][candidate]
		#	#print path[-10:]
		if cost==unused_cost:
			if unused_cost>max_so_far[0]+100 or unused_cost>-100:
				print unused_cost#,path
				max_so_far[0]=unused_cost
			local_path=deepcopy(path)
			local_path.append(candidate)
			results.append((len(local_path),flr,local_path,cost))
		else:
			print >> sys.stderr, "ERROR !!!" , cost
			sys.exit(1)
	return results,False
	#for result in results:
	#	local_flr=deepcopy(fl)
	#	local_flr[last_node][result]-=1
	#	#print flr,path
	#	local_path=deepcopy(path)
	#	local_path.append(result)
	#	search(local_flr,local_path)
	

params['walks']=read_graph(sys.stdin)
fname='/dev/shm/out_%d' % pid
print print_graph(fname)
flx,cx=get_flow(fname)
flx_keys=flx.keys()
random.shuffle(flx_keys)
start_node=flx_keys[0]

#tflow=setup()

#p = Pool(32)

def look(start,end):
	print_graph(fname)
	flx,cost=get_flow(fname)
	q=[(1,{},[start_node],cost)]
	while len(q)>0:
		tpl=q.pop()
		results,d = search(tpl)
		if d or time.time()>end:
			return
		for result in results:
			q.append(result)


origin_time = time.time()
max_time=60*60*60

d=60*60
while True:
	start = time.time()
	look(start,start+d)
	if time.time()<start+d:
		#found a solution
		print "FOUND A SOLUTION"
	else:
		#didnt find one
		#d=int(d*1.2)
		print "NO SOLUTION FOUND"
	print >> sys.stderr, "Restarting with ", d
