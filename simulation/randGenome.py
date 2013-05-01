#!/usr/bin/python
import random
import sys



alphabet=['A','C','G','T']

def usage():
	print sys.argv[0], " length"


def get(length,r):
	ret=""
	while length>0:
		ret+=r.choice(alphabet)
		length-=1
	return ret

if __name__=='__main__':
	if len(sys.argv)!=2:
		usage()
		sys.exit(1)
	length=int(sys.argv[1])
	r=random.Random()
	
	perline=79

	print ">RANDOM_GENOME_%d" % length
	while length>0:
		if length>=perline:
			print get(perline,r)
			length-=perline
		else:
			print get(length,r)
			length=0
	
	
