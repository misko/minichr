

import sys
import gzip


if len(sys.argv)!=3:
	print "%s problemxfile.gz edge_lookup" % sys.argv[0]
	sys.exit(1)

pb_fname=sys.argv[1]
el_fname=sys.argv[2]

nodes={}
ins={}
outs={}
edges=[]


'''
Somatic edges types
	0 + + >->
	1 - - <-<
	2 + - >-<
	3 - + <->
'''

edge_lookup=[]

def addxedge(fn,tn,ty,cost,cap):
	eid=len(edges)
	edge_lookup.append((eid,fn,tn))
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

def readxpb(fname):
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
				addxedge(fn,tn,0,cost,cap)
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
				addxedge(fn,tn,ty,cost,cap)
			elif line[2]=='N':
				#its a node
				#c NODE 19 chr1:808212
				line=line.split()
				nid=int(line[2])
				nodes[line[3]]=nid
				ins[nid]=[]
				outs[nid]=[]
	h.close()




readxpb(pb_fname)

def outputxlpxsolve():
	objxfunc=[]
	for x in range(len(edges[e])):
			if edges[e][x][0]!=0:
				objxfunc.append("%d %s" % (edges[x][0],eid))
	print "min: "," + ".join(objxfunc),";"
	for n in nodes:
		c=[]
		for x in ins[nodes[n]]:
			c.append(" + f%d" % (x))
		for x in outs[nodes[n]]:
			c.append(" - f%d" % (x)) 
		if len(c)!=0:
			print "".join(c),"= 0",";"
	for e in edges:
		for x in range(len(edges[e])):
			eid="%sx%d" % (e,x)
			print " %s <= %d" % (eid,edges[eid][x][1]),';'
	print "int "
	o=[]
	for e in edges:
		for x in range(len(edges[e])):
			eid="%sx%d" % (e,x)
			o.append(eid)
	print ",".join(o),';'

def outputxcplex():
	#now output the lp file
	objxfunc=[]
	for x in range(len(edges)):
		if edges[x][0]!=0:
			objxfunc.append("%d f%d" % (edges[x][0],x))
	print "Minimize\n obj:",
	ll=5
	for x in objxfunc:
		if ll+len(x)+3>=80:
			ll=0
			print ""
			print " +",x,
			ll+=len(x)+3
		else:
			print " +",x,
			ll+=len(x)+3

	print "Subject To"
	for n in nodes:
		c=[]
		for x in ins[nodes[n]]:
			c.append(" + f%d" % (x))
		for x in outs[nodes[n]]:
			c.append(" - f%d" % (x)) 
		if len(c)!=0:
			print " n%d: " % nodes[n],"".join(c)," = 0"
			#print "".join(c),"= 0",";"

	print "Bounds"
	for x in range(len(edges)):
		print " f%d <= %d" % (x,edges[x][1])
	print "Integers "
	o=[]
	for x in range(len(edges)):
		o.append("f%d" % x)
	print " ".join(o)
	print "End"

outputxcplex()

h=open(el_fname,'w')
h.write("\n".join(map(lambda x : str(x) , edge_lookup)))
h.close()
