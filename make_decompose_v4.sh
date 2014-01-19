#!/bin/bash

g=/filer/misko/mini_chr/git/minichr/
p=/filer/misko/picard/picard-tools-1.56/
s=/home/misko/apps/bin/samtools
j=/usr/bin/java
b=/filer/misko/bedtools-2.17.0/
c=/data/misko/2013.04.12/cs2-4.3/cs2.exe
export LD_LIBRARY_PATH=/home/buske/arch/sge6/lib:/home/buske/arch/sge6/lib::/dupa-filer/misko/gurobi/gurobi550/linux64/lib/:/dupa-filer/misko/gurobi/gurobi550/linux64/lib/
#using hg19
#ref=/filer/hg19/hg19.fa
#or hg18
#ref=/filer/hg18/hg18.fa

if [ $# -ne 1 ]; then 
	echo $0 folder
	exit
fi
wd=$1

function decomp {
	
	mkdir -p $wd/decomps
	pushd $wd/decomps
	dname=LPY`tr -cd '[:alnum:]' < /dev/urandom | fold -w30 | head -n1`
	zcat ../Qproblem_file_q1350sq300_m3.solved.gz > ${dname}.tmp
	/usr/bin/pypy $g/decompose_flow_v2.py ${dname}.tmp $dname
	mv $dname ${dname}.lp
	rm ${dname}.tmp
	popd 
}


x=0
while [ $x -lt 10 ]; 
do
	decomp
	x=`expr $x + 1`
done
