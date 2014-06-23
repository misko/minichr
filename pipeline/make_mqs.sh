#!/bin/bash

g=/filer/misko/mini_chr/git/minichr/
p=/filer/misko/picard/picard-tools-1.56/
s=/home/misko/apps/bin/samtools
j=/usr/bin/java
b=/filer/misko/bedtools-2.17.0/
c=/data/misko/2013.04.12/cs2-4.3/cs2.exe

#using hg19
#ref=/filer/hg19/hg19.fa
#or hg18
#ref=/filer/hg18/hg18.fa

if [ $# -ne 1 ]; then 
	echo $0 folder
	exit
fi
wd=$1

function moverlap {
	pushd $wd
	ref=`python $g/get_ref.py downloadable.xml | awk '{print $NF}' | head -n 1`
	$g/mapability/bigWigAverageOverBed $g/mapability/$ref/* nsubtract_centrosubtract_1000bp.bed nsubtract_centrosubtract_1000bp.bed.out
	cat nsubtract_centrosubtract_1000bp.bed.out | awk '{print $1,$NF}' > nsubtract_centrosubtract_1000bp_mqs 
	popd
}

#run overlaps
moverlap




