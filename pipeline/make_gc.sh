#!/bin/bash

g=/filer/misko/mini_chr/git/minichr/
p=/filer/misko/picard/picard-tools-1.56/
s=/home/misko/apps/bin/samtools
j=/usr/bin/java
b=/filer/misko/bedtools-2.17.0/
c=/data/misko/2013.04.12/cs2-4.3/cs2.exe


if [ $# -ne 1 ]; then 
	echo $0 folder
	exit
fi
wd=$1

function moverlap {
	pushd $wd

	ref=`python $g/get_ref.py downloadable.xml | awk '{print $NF}' | head -n 1 | awk '{print $NF}'`
	if [ -z "$ref" ] ;  then
		echo "FAILED!!! COULD NOT FIND REF"
		touch fail
		exit
	fi
	echo $ref > ref
	$g/getcov/print_cov ${wd}/tumor_cov.gz /dupa-filer/${ref}.fa.gz > ${wd}/tumor_cov.gz.log
	$g/getcov/print_cov ${wd}/normal_cov.gz /dupa-filer/${ref}.fa.gz > ${wd}/normal_cov.gz.log
	# | sort | uniq > nsubtract_centrosubtract_1000bp_mqs
	popd
}

#run overlaps
moverlap




