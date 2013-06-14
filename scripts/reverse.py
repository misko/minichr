#!/usr/bin/python


import sys

rev={'N':'N','n':'n','A':'T','a':'t','C':'G','c':'g','G':'C','g':'c','T':'A','t':'a'}

def reverse(s):
	ret=""
	s=s.strip()
	for x in range(len(s)):
		ret+=rev[s[-x-1]]
	return ret


if __name__=='__main__':
	if len(sys.argv)==1:
		print '%s [-f filename] | string' % sys.argv[0]
		sys.exit(1)
	if len(sys.argv)==2:
		print reverse(sys.argv[1])
	else:
		h=open(sys.argv[2])
		lines=h.readlines()
		prev=0
		title=""
		for x in range(len(lines)):
			cur=lines[x]
			if cur[0]=='>':
				#got a title start recording, dump previous
				if prev!=0:
					#dump previous
					subset=map(reverse,lines[prev:x])
					subset.reverse()
					print title.strip()
					print "\n".join(subset)
				prev=x+1
				title=cur
		#check if title non-empty
		if len(title)!=0:
			subset=map(reverse,lines[prev:x+1])
			subset.reverse()
			print title.strip()
			print "\n".join(subset)
				

