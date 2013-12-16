#!/bin/bash

export PATH=$PATH:/data/misko/2013.04.12/cs2-4.3/


ym=1
y=1
while [ $y -le $ym ]; do
	echo "running with new flow " $y
	./walker /dupa-filer/misko/dipg29t/n3t10_nocentro_refined_new /dupa-filer/misko/dipg29t/n3t10_nocentro_refined_breakpoints_normal_mq10 ../out15 $y


	mkdir outputs_$y	
	cp problem_file outputs_$y	

	xm=32
	x=1
	while [ $x -le $xm ]; do
	 	cat problem_file | pypy path_finder.py  > outputs_$y/out_solved_$x  &
		x=`expr $x + 1`
	done


	wait

	y=`expr $y + 1`
done


#xs="10 11 12 13 14 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42"
#for x in $xs; do
#	 cat problem_file | pypy path_finder.py  > outputs/out_solved_$x  &
#	 sleep 1
#done
