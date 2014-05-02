import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import sys


def get_fractions(f):
	h=open(f)
	counts=[]
	lines=h.readlines()
	for x in xrange(len(lines)):
		line=lines[x].strip()
		if line=="Finding the average" and x<len(lines)-1:
			try:
				counts=map(lambda x : int(x), lines[x+1].strip().split())
				s=sum(counts)
				for x in xrange(len(counts)):
					counts[x]/=float(s)
			except:
				pass
	h.close()
	return counts

if len(sys.argv)!=4:
	print "%s normal tumor fig_file" % sys.argv[0]
	sys.exit(1)

normal_filename=sys.argv[1]
tumor_filename=sys.argv[2]
fig_filename=sys.argv[3]

if fig_filename[-3:]!='png':
	print >> sys.stderr, "MUST BE png filetype!"
	sys.exit(1)


normal=get_fractions(normal_filename)
tumor=get_fractions(tumor_filename)

diffs=[]

for x in range(len(normal)):
	diffs.append(abs(normal[x]-tumor[x]))


x = range(len(diffs))
plt.scatter(x, diffs,color='black')
plt.scatter(x, tumor,color='red')
plt.scatter(x, normal,color='blue')
plt.savefig(fig_filename)

print sum(diffs)
