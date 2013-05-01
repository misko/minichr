#!/bin/bash


reads=700

runs="type0 type1 type1_small type2_3"

for run in $runs; do
	echo $run
	python simulate_reads.py $run/don.fa $reads
	mv ls_?_reads.fq $run/
	pushd $run
	sh commands
	popd
	samtools view $run/mapped_sorted.bam | ../clustering/cluster 700 70 > $run/links
	cat $run/links | sed 's/\([0-9]\)	\([0-9]*\)[:]\([0-9]*\)	\([0-9]*\)[:]\([0-9]*\)	\([0-9]*\)	\(.*\)/\2	\3	\4	\5	\1	\6/g' | awk '{if ($1==23) {$1="X"}; if ($1==24) {$1="Y"}; if ($3==23) {$3="X"}; if ($3==24) {$3="Y"}; a="+"; b="+"; if ($5==1) {a="-"; b="-"}; if ($5==2) {a="+"; b="-"}; if ($5==3) {a="-"; b="+"}; print "chr"$1"\t"$2"\t"a"\tchr"$3"\t"$4"\t"b"\t"$6}' > $run/links_f
	
done

