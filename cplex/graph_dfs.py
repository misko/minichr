

import sys
import gzip


if len(sys.argv)!=3:
	print "%s problemxfile.gz lp_solution" % sys.argv[0]
	sys.exit(1)

pb_fname=sys.argv[1]
lp_sol_fname=sys.argv[2]

nodes={}
ins={}
outs={}
edges=[]
contigs=[0]

'''
Somatic edges types
	0 + + >->
	1 - - <-<
	2 + - >-<
	3 - + <->
'''

edge_lookup=[]

def add_edge(fn,tn,ty,cost,cap):
	eid=len(edges)
	edge_lookup.append((eid,fn,tn,ty))
	edges.append((cost,cap))
	if ty==0:
		outs[fn].append(eid)
		ins[tn].append(eid)
	elif ty==1:
		ins[fn].append(eid)
		outs[tn].append(eid)
	elif ty==2:
		outs[fn].append(eid)
		outs[tn].append(eid)
	elif ty==3:
		ins[fn].append(eid)
		ins[tn].append(eid)


def read_pb(fname):
	h=gzip.open(fname)
	for line in h:
		if line[0]=='c':
			if line[2]=='G':
				#c Genomic       chr1:712095     chr1:714492     66      4       2397    736     1844
				#ss << "c Genomic\t" << e.posa.str() << "\t" << e.posb.str()  << "\t" << cost << "\t" << cap << "\t" << e.length() << "\t" << ei.normal << "\t" << ei.tumor << endl;
				line=line.split()
				cost=int(line[4])
				cap=int(line[5])
				l=int(line[6])
				normal=int(line[7])
				tumor=int(line[8])
				fn=nodes[line[2]]
				tn=nodes[line[3]]
				add_edge(fn,tn,0,cost,cap)
			elif line[2]=='S':
				#c Somatic       chr1:814892     chrY:11934477   1       23      399     237     33
				#ss << "c Somatic\t" << e.posa.str() << "\t" << e.posb.str() << "\t" << ei.type << "\t" << cost << "\t" << cap << "\t" << ei.normal << "\t" << ei.tumor << endl;
				line=line.split()
				ty=int(line[4])
				cost=int(line[5])
				cap=int(line[6])
				normal=int(line[7])
				tumor=int(line[8])
				fn=nodes[line[2]]
				tn=nodes[line[3]]
				#if line[2]=="chrX:154507287" and line[3]=="chrX:154510610":
				#	print line,fn,tn
				#	sys.exit(1)
				#if fn==257490 and tn==257493:
				#	print line
				#	print len(edges)
				#	sys.exit(1)
				add_edge(fn,tn,ty,cost,cap)
			elif line[2]=='N':
				#its a node
				#c NODE 19 chr1:808212
				line=line.split()
				nid=int(line[2])
				if nid%2==0:
					nodes[line[3]]=nid
					ins[nid]=[]
					outs[nid]=[]
		elif line[0]=='a' and line[2]=='1' and line[3]=='\t' and line[4]=='3':
			print >> sys.stderr, line
			contigs[0]=int(line[6])
			print >> sys.stderr, "Considering ", contigs[0], "contigs"
	#add in the funky edges	
	outs[0]=[]
	ins[0]=[]
	ins[2]=[]
	outs[2]=[]
	add_edge(2,0,0,-10,1000)
	for node in nodes:
		nid=nodes[node]
		add_edge(0,nid,0,5,1000)
		add_edge(0,nid,2,5,1000)
		add_edge(nid,2,0,5,1000)
		add_edge(nid,2,3,5,1000)
	h.close()




read_pb(pb_fname)


def dfs(red_edges, found):
	rnodes=

