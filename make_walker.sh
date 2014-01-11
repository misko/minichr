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

function mwalker {
	i=$1
	sq=$2
	m=$3
	pushd $wd
	$g/walker/walker_c $wd/nsubtract_centrosubtract_1000bp ${wd}/hmm N 0 ${m} 2 > walker_out_q${i}sq${sq}_m${m}
	mv problem_file.gz problem_file_q${i}sq${sq}_m${m}.gz
	zcat  problem_file_q${i}sq${sq}_m${m}.gz | /data/misko/2013.04.12/cs2-4.3/cs2.exe | gzip > solved_q${i}sq${sq}_m${m}.gz

	#make the lp problem
	pypy /filer/misko/mini_chr/git/minichr/cplex/graph_c_to_lp.py problem_file_q${i}sq${sq}_m${m}.gz el_q${i}sq${sq}_m${m}.gz nsubtract_centrosubtract_1000bp_mqs ${i} ${sq} ${m} > lp_prob_q${i}sq${sq}_m${m}.lp
	#sovle LP
	#/dupa-filer/misko/gurobi/gurobi550/linux64/bin/gurobi_cl MIPGap=0.15 Threads=8 ResultFile=result_q${i}sq${sq}_m${m}.sol lp_prob_q${i}sq${sq}_m${m}.lp
	/dupa-filer/misko/gurobi/gurobi550/linux64/bin/gurobi_cl MIPGap=0.15 Threads=8 ResultFile=result_q${i}sq${sq}_m${m}.sol lp_prob_q${i}sq${sq}_m${m}.lp
	#conver LP to flow
	pypy /filer/misko/mini_chr/git/minichr/cplex/lp_to_flow.py el_q${i}sq${sq}_m${m}.gz result_q${i}sq${sq}_m${m}.sol | gzip > lp_flow_q${i}sq${sq}_m${m}.gz
	#flow to graph
	pypy /filer/misko/mini_chr/git/minichr/walker/flow_to_graph_v2.py problem_file_q${i}sq${sq}_m${m}.gz lp_flow_q${i}sq${sq}_m${m}.gz ${m} 0 1000 > g_lp_q${i}sq${sq}_m${m}
	#update graph with no flows
	pypy /filer/misko/mini_chr/git/minichr/insert_no_flows.py g_lp_q${i}sq${sq}_m${m} hmm > g_lp_q${i}sq${sq}_m${m}_whmm
	popd 
}

#run overlaps
#mwalker 4000 1000
mwalker 1350 300 3 # contigQ somatic_base_Q multiplier 

