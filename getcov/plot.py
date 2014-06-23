#!/usr/bin/python
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import sys


def get_fractions(f):
	h=open(f)
	gc_counts=[]	
	rpc_counts=[]
	cov=0
	for l in h:
		l=l.strip().split()
		if l[0]=="GC":
			gc_counts=map(lambda x : int(x), l[1:])
			s=sum(gc_counts)
			for x in xrange(len(gc_counts)):
				if s!=0:
					gc_counts[x]/=float(s)
				else:
					gc_counts[0]=0
		elif l[0]=="RPC":
			rpc_counts=map(lambda x : int(x), l[1:-1])
			s=sum(rpc_counts)
			cov=s
			for x in xrange(len(rpc_counts)):
				if s!=0:
					rpc_counts[x]/=float(s)
				else:
					rpc_counts[x]=0
	h.close()
	return (gc_counts,rpc_counts,cov)

if len(sys.argv)!=4:
	print "%s normal tumor fig_file" % sys.argv[0]
	sys.exit(1)

normal_filename=sys.argv[1]
tumor_filename=sys.argv[2]
fig_filename=sys.argv[3]

if fig_filename[-3:]!='png':
	print >> sys.stderr, "MUST BE png filetype!"
	sys.exit(1)


gc_normal,rpc_normal,cov_normal=get_fractions(normal_filename)
gc_tumor,rpc_tumor,cov_tumor=get_fractions(tumor_filename)

diffs=[]

for x in range(len(gc_normal)):
	diffs.append(abs(gc_normal[x]-gc_tumor[x]))


x = range(len(diffs))
plt.scatter(x, diffs,color='black')
plt.scatter(x, gc_tumor,color='red')
plt.scatter(x, gc_normal,color='blue')
plt.savefig(fig_filename)
print sum(diffs),min(rpc_normal),min(rpc_tumor),cov_normal,cov_tumor,max(float(cov_tumor)/float(cov_normal),float(cov_normal)/float(cov_tumor))
