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

if [ $# -ne 2 ]; then 
	echo $0 folder tcga-id
	exit
fi
wd=$1
id=$2

function mwalker {
	i=$1
	sq=$2
	m=$3
	pushd $wd	

	ref=`python $g/get_ref.py  downloadable.xml  | awk '{print $NF}' | head -n 1`
	group=`$g/get_group.sh downloadable.xml`
	echo $id $group > group
	/usr/bin/pypy $g/walker/contigs_to_coords.py Qproblem_file_q${i}sq${sq}_m${m}.gz decompositions_q${i}sq${sq}_m${m}.loops.full ip_prob_q${i}sq${sq}_m${m}.lp.sol ${m} $g/genes/genes_wstrands_formatted_${ref}.txt $g/genes/oncos ${id} group > stats

	popd 
}

#run overlaps
#mwalker 4000 1000

#mwalker 1350 300 3 # contigQ somatic_base_Q multiplier 

#this was run with ip SQ = SQ / 2
#mwalker 1350 290 3 # contigQ somatic_base_Q multiplier 

#this was run with ip SQ = -20
#mwalker 1350 291 3 # contigQ somatic_base_Q multiplier 

#this was run with ip SQ = 20
mwalker 1350 292 3 # contigQ somatic_base_Q multiplier 
