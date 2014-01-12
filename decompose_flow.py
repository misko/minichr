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
count=0

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
		count+=fl
	

for x in ins:
	if ins[x]!=outs[x]:
		print >> sys.stderr, "FAIL",x
		sys.exit(1)

for x in outs:
	if ins[x]!=outs[x]:
		print >> sys.stderr, "FAIL",x
		sys.exit(1)
		

def get_start_edge():
	es=[]
	for fn in edges:
		for tn in edges[fn]:
			if edges[fn][tn]<=0:
				print >> sys.stderr, "FAILDE!"
				sys.exit(1)
			for x in range(edges[fn][tn]):
				es.append((fn,tn))
	if len(es)==0:
		return False
	return choice(es)			

def decrease(fn,tn):
	edges[fn][tn]-=1
	if edges[fn][tn]==0:
		edges[fn].pop(tn,None)
	if len(edges[fn])==0:
		edges.pop(fn,None)



#find the neighbours
stack=[]
while True:
	if len(stack)<2:
		start_edge=get_start_edge()
		if start_edge==False:
			break
		decrease(start_edge[0],start_edge[1])
		count-=1
		stack=[start_edge[0],start_edge[1]]
	fn=stack[-1]
	ns=[]
	for tn in edges[fn]:
		for x in range(edges[fn][tn]):
			ns.append(tn)
	if len(ns)==0:
		break
	#choose a neighbour 
	tn=choice(ns)
	decrease(fn,tn)
	count-=1
	if tn in stack:
		#remove the loop!
		idx=stack.index(tn)
		loop=stack[idx:]
		print loop
		stack=stack[:idx]
	stack.append(tn)



	
