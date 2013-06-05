#!/usr/bin/bash

g=/filer/misko/mini_chr/git/minichr/
p=/filer/misko/picard/picard-tools-1.56/
s=samtools
j=java
b=/filer/misko/bedtools-2.17.0/
ref=/filer/hg19/hg19.fa


if [ $# -ne 3 ]; then 
	echo $0 normal_bam tumor_bam work_directory
	exit
fi

normalbam=$1
if [ ! -f $normalbam ]; then
	echo Cannot find normal bam file!
	exit
fi
if [ ! -f $normalbam.bai ]; then
	echo Bam file not indexed!
	exit
fi
tumorbam=$2
if [ ! -f $tumorbam ]; then
	echo Cannot find normal bam file!
	exit
fi
if [ ! -f $tumorbam.bai ]; then
	echo Bam file not indexed!
	exit
fi
wd=$3
if [ -d ${wd} ] ; then
	echo Temp directory exists!
	exit
fi

mkdir -p $wd

#link in the files
ln -s $normalbam $wd/normal.bam
ln -s ${normalbam}.bai $wd/normal.bam.bai
ln -s $tumorbam $wd/tumor.bam
ln -s ${tumorbam}.bai $wd/tumor.bam.bai


tys="normal tumor"
cluster=$g/clustering/cluster
cov_mapq=20
cluster_mapq=10
for ty in $tys; do
	#BAM file coverage
	echo Skipping cov file generation
	$s mpileup -q ${cov_mapq} $wd/${ty}.bam | $g/getcov/get_cov ${wd}/${ty}_cov
	echo "Warning sampling first 100000 for insert size"
	$j -jar $p/CollectInsertSizeMetrics.jar I=$wd/${ty}.bam H=$wd/${ty}_histo STOP_AFTER=100000 O=$wd/${ty}_stats
	cat $wd/${ty}_stats | grep -A 1 MED | grep -v "\-\-" | grep -v MEDIAN | awk '{m+=$5; s+=$6} END {print m/NR,s/NR}' > $wd/${ty}_mean_and_std.txt
	mean=`cat ${wd}/${ty}_mean_and_std.txt | awk '{print int($1)}'`
	stddev=`cat ${wd}/${ty}_mean_and_std.txt | awk '{print int($2)}'`
	
	echo Skipping chrM 
	$s view -q ${cluster_mapq} ${wd}/${ty}.bam | grep -v chrM | head -n 100000 | $cluster $mean $stddev | gzip > ${wd}/${ty}_clusters.txt.gz
	covs="10 5 3 2 0"
	for cov in $covs; do
		zcat ${wd}/${ty}_clusters.txt.gz | awk '{if ($4>10) {print $0}}' | sed 's/\([0-9]\)\t\([0-9]*\)[:]\([0-9]*\)\t\([0-9]*\)[:]\([0-9]*\)\t\([0-9]*\)\t\(.*\)/\2\t\3\t\4\t\5\t\1\t\6/g' | awk '{if ($1==23) {$1="X"}; if ($1==24) {$1="Y"}; if ($3==23) {$3="X"}; if ($3==24) {$3="Y"}; a="+"; b="+"; if ($5==1) {a="-"; b="-"}; if ($5==2) {a="+"; b="-"}; if ($5==3) {a="-"; b="+"}; print "chr"$1"\t"$2"\t"a"\tchr"$3"\t"$4"\t"b"\t"$6}' |  gzip > ${wd}/${ty}_cov${cov}.txt_f.gz
	done

done

python2.6 $g/scripts/overlap.py 1000 ${wd}/normal_cov0.txt_f.gz ${wd}/tumor_cov5.txt_f.gz  |  awk '{if (NF>4) {print $0}}' | awk '{print $1"\t"$2"\t"$2+1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}'  | $b/bin/intersectBed -a - -b $g/centromeres_merged -v | awk '{print $5"\t"$6"\t"$6+1"\t"$7"\t"$1"\t"$2"\t"$4"\t"$8}' | $b/bin/intersectBed -a - -b $g/centromeres_merged -v  | awk '{print $5"\t"$6"\t"$7"\t"$1"\t"$2"\t"$4"\t"$8}' | awk '{d=$2-$5; if (d<0) {d=-d}; if ($3==$6 && $1==$4 && d<2000) { } else {print $0}}' > $wd/nsubtract_centrosubtract_1000bp




