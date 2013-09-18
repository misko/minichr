

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
contigs=[0]
edge_by_nodes={}

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
	if fn<0 or tn<0:
		print >> sys.stderr, "Failed to use nodes properly"
		sys.exit(1)
	k=(fn,tn,ty)
	if k not in edge_by_nodes:
		edge_by_nodes[k]={}
	if cost not in edge_by_nodes[k]:
		edge_by_nodes[k][cost]=[]
	edge_by_nodes[k][cost].append(eid)
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
	#balance
	for n in nodes:
		c=[]
		for x in ins[nodes[n]]:
			c.append(" + f%d" % (x))
		for x in outs[nodes[n]]:
			c.append(" - f%d" % (x)) 
		if len(c)!=0:
			print " n%d: " % nodes[n],"".join(c)," = 0"
			#print "".join(c),"= 0",";"
	#put in the source and sink
	cx=[]
	for w in range(contigs[0]):
		cx.append(" - 10") # - sum(m_i) TODO
	c=[]
	for x in outs[0]:
		c.append(" + f%d" % (x)) 
	if len(c)!=0:
		print " n%d: " % nodes[n],"".join(c),"".join(cx)," =",0
	c=[]
	for x in ins[0]:
		c.append(" + f%d" % (x)) 
	if len(c)!=0:
		print " n%d: " % nodes[n],"".join(c),"".join(cx)," =",0
	c=[]
	for x in outs[2]:
		c.append(" + f%d" % (x)) 
	if len(c)!=0:
		print " n%d: " % nodes[n],"".join(c),"".join(cx)," =",0
	c=[]
	for x in ins[2]:
		c.append(" + f%d" % (x)) 
	if len(c)!=0:
		print " n%d: " % nodes[n],"".join(c),"".join(cx)," =",0
	#also the meta edges
	mid={} # fn,tn,ty meta edge id
	for fn,tn,ty in edge_by_nodes:
		c=[]
		mid[(fn,tn,ty)]=len(mid)
		for w in range(contigs[0]):
			c.append(" + 10 w%dx%d" % (w,mid[(fn,tn,ty)])) # TODO m_i * w_i_eid
		for cost in edge_by_nodes[(fn,tn,ty)]:
			for eid in edge_by_nodes[(fn,tn,ty)][cost]: 
				if eid<0:
					print >> sys.stderr, "ERROR"
					sys.exit(1)
				c.append(" - f%d" % eid)
		print " w%d: " % mid[(fn,tn,ty)],"".join(c)," =",0
	
	
	print "Bounds"
	for x in range(len(edges)):
		print " f%d <= %d" % (x,edges[x][1])
	print "Integers "
	o=[]
	for x in range(len(edges)):
		o.append("f%d" % x)
	#for x in range(contigs[0]):
	#	o.append("m%d" % x)
	for fn,tn,ty in edge_by_nodes:
		for w in range(contigs[0]):
			o.append("w%dx%d" % ( w,mid[(fn,tn,ty)]))
	print " ".join(o)
	print "End"

outputxcplex()

#check for straight aways
sins={}
souts={}
for fn,tn,ty in edge_by_nodes:
	if fn not in souts:
		souts[fn]=set()
	souts[fn].add(tn)
	if tn not in sins:
		sins[tn]=set()
	sins[tn].add(fn)
al=set()
al.update(set(sins.keys()))
al.update(set(souts.keys()))
for n in al:
	if n in souts and n in sins:
		if len(souts[n])==1 and len(sins[n])==1:
			print >> sys.stderr,  "Straight through edge ", n 

#write the output
h=gzip.open(el_fname,'w')
h.write("\n".join(map(lambda x : str(x) , edge_lookup)))
h.close()
