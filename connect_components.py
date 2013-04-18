#!/usr/bin/python


import sys


if len(sys.argv)!=3:
	print "%s links segmentation" % sys.argv[0]
	sys.exit(1)



links_filename=sys.argv[1]
seg_filename=sys.argv[2]

genomic_edges={}
cluster_edges={}
arc_edges={}
nodes=[]
pos_to_node={}
starts={}
ends={}

def to_chr(s):
	s=s.lower()
	if s[:3]=='chr':
		s=s[3:]
	if s=='x':
		return 23
	if s=='y':
		return 24
	if s=='m':
		return 25
	return int(s)


#read in segments
last_chr=0
seg_file=open(seg_filename,'r')
for line in seg_file:
	if line[0]=="#":
		continue
	#cp	start	end	length	cancer	lambda(normal/2)
	#3       chr1:60457326   chr1:60459726   2400    846     318
	line=line.split()
	cp=int(line[0])
	posa=(to_chr(line[1].split(':')[0]),int(line[1].split(':')[1]))
	if posa not in nodes:
		pos_to_node[posa]=len(nodes)
		nodes.append(posa)
	posb=(to_chr(line[2].split(':')[0]),int(line[2].split(':')[1]))
	if posb not in nodes:
		pos_to_node[posb]=len(nodes)
		nodes.append(posb)
	cancer=int(line[4])
	lmbda=int(line[5])
	starts[posa]=(posa,posb)
	ends[posb]=(posa,posb)
	genomic_edges[(posa,posb)]={'cp':cp,'posa':posa,'posb':posb,'cancer':cancer,'lmbda':lmbda}

#read in the links
links_file=open(links_filename,'r')
for line in links_file:
	if line[:3]!='chr':
		continue
	#chr1    201862828       205356414       1       124     0       0       0       0       chr1    EDGE
	line=line.split()
	posa=(to_chr(line[0]),int(line[1]))
	if posa not in pos_to_node:
		print posa
		print "Fatal error posa"
		sys.exit(1)
	posb=(to_chr(line[9]),int(line[2]))
	print posb
	if posb not in pos_to_node:
		print "Fatal error posb"
		sys.exit(1)
	t=int(line[3])
	if t%2==0:
		#goes in on positve
		#check if previous edge has cp>=3
		e=ends[posa]
		if genomic_edges[e]['cp']<3:
			continue
		if t==0:
			#leaves on positive
			e=starts[posb]
			if genomic_edges[e]['cp']<3:
				continue
		else:
			#leaves on negative
			e=ends[posb]
			if genomic_edges[e]['cp']<3:
				continue
	else:
		#goes in on negativep
		e=starts[posa]
		if genomic_edges[e]['cp']<3:
			continue
		if t==3:
			#leaves on positive
			e=starts[posb]
			if genomic_edges[e]['cp']<3:
				continue
		else:
			#leaves on negative		
			e=ends[posb]
			if genomic_edges[e]['cp']<3:
				continue
	cluster_edges[(posa,posb)]={'type':int(line[3]),'support':int(line[4])}



import networkx as nx
G=nx.Graph()

for posa,posb in genomic_edges:
	d=genomic_edges[(posa,posb)]
	if d['cp']>=3:
		G.add_node(posa)
		G.add_node(posb)
		G.add_edge(posa,posb)

for posa,posb in cluster_edges:
	G.add_edge(posa,posb)

components=nx.connected_components(G)
s_components=[]
for component in components:
	if len(component)>2:
		s_components.append((len(component),component))
s_components.sort(reverse=True)

print s_components[-3:]

nodes=set()
for component in s_components[:1]:
	for node in component[1]:
		nodes.add(node)


Gx=nx.Graph()
xgenomic_edges=[]
for posa,posb in genomic_edges:
	d=genomic_edges[(posa,posb)]
	if d['cp']>=3 and (posa in nodes) and (posb in nodes):
		Gx.add_node(posa)
		Gx.add_node(posb)
		Gx.add_edge(posa,posb,weight='1')
		xgenomic_edges.append((posa,posb))



import matplotlib.pyplot as plt
pos=nx.spring_layout(Gx) # positions for all nodes

print pos


xcluster_edges=[]
for posa,posb in cluster_edges:
	if (posa in nodes) and (posb in nodes):
		Gx.add_edge(posa,posb,weight='1')
		xcluster_edges.append((posa,posb))

# nodes
print pos
nx.draw_networkx_nodes(Gx,pos)
labels={}
for node in nodes:
	labels[node]=str(node)
	nx.draw_networkx_labels(G,pos,labels,font_size=12)

for posa,posb in xgenomic_edges:
	cp=5+genomic_edges[(posa,posb)]['cancer']/(1+genomic_edges[(posa,posb)]['lmbda'])
	nx.draw_networkx_edges(Gx,pos,
	                       edgelist=[(posa,posb)],
	                       width=cp/2,alpha=0.5,edge_color='r')
for posa,posb in xcluster_edges:
	nx.draw_networkx_edges(Gx,pos,
	                       edgelist=[(posa,posb)],
	                       width=cluster_edges[(posa,posb)]['support']/2+2.5,alpha=0.5,edge_color='b')
	
# edges
#nx.draw_networkx_edges(G,pos,width=1.0,alpha=0.5)
#nx.draw_networkx_edges(G,pos,
#                       edgelist=[(0,1),(1,2),(2,3),(3,0)],
#                       width=8,alpha=0.5,edge_color='r')
#nx.draw_networkx_edges(G,pos,
#                       edgelist=[(4,5),(5,6),(6,7),(7,4)],
#                       width=8,alpha=0.5,edge_color='b')



#nx.draw_spring(Gx)
#plt.show()
plt.savefig('graph.png')


