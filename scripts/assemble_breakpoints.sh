#!/bin/bash

s=samtools
g=/filer/misko/mini_chr/git/minichr/
p=/filer/misko/picard/picard-tools-1.56/
v=/home/misko/apps/bin/
n=/home/misko/naive/naive
c=/home/misko/apps/bin/chaos
b=blat

if [ $# -ne 6 ]; then
	echo $0 ref links bam mean stddev temp
	exit
fi


r=$1
links=$2
bam=$3
mean=$4
stddev=$5
wd=$6

if [ ! -f $links ] ; then
	echo Links file $links does not exist!
	exit
fi

if [ ! -f $bam ] ; then
	echo Bam file $bam does not exist!
	exit
fi

if [ -d $wd ] ; then
	echo temp directory exists!
	#exit
fi

mkdir -p $wd


cat $links | awk '{OFS="\t"; type=0; if ($3=="+") {if ($6=="+") {type=0} else {type=2} } else { if ($6=="+") {type=3} else {type=1} }; if ($1>$4 || ($1==$4 && $2>$5) ) { if (type<2) {type=1-type} ; fc=$4; fcc=$5; tc=$1; tcc=$2;} else {fc=$1; fcc=$2; tc=$4; tcc=$5;} ; print type,fc,fcc,tc,tcc}' | sort | uniq > $wd/uniq_links


#cat $links | awk '{if ($1>$4 || ($1==$4 && $2>$5)) {if ($3=="+") {$3="-"} else {$3="+"}; if ($6=="+") {$6="-"} else {$6="+"}; print $4"\t"$5"\t"$6"\t"$1"\t"$2"\t"$3} else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}'  | sort | uniq > $wd/uniq_links

#cat $wd/uniq_links | awk -v f=1000 -v b=300 '{if ($3=="-") {as=$2-b; ae=$2+f} else {as=$2-f; ae=$2+b}; if ($6=="-") {bs=$5-f; be=$5+b} else {bs=$5-b; be=$5+f} ; if (as<=0) {as=1}; if (bs<=0) {bs=1}; print $1":"as"-"ae"\t"$4":"bs"-"be}' > $wd/regions


$s view -H $bam > $wd/sam_header

i=0

while read rawline; do 
	echo working on $rawline
	mkdir -p $wd/$i
	echo $rawline > $wd/$i/line
	line=`echo $rawline | awk -v f=5000 -v b=5000 '{if ($1==0 || $1==2) {as=$3-f; ae=$3+b} else {as=$3-b; ae=$3+f}; if ($1==0 || $1==3) {bs=$5-b; be=$5+f} else {bs=$5-f; be=$5+b} ; if (as<=0) {as=1}; if (bs<=0) {bs=1}; print $2":"as"-"ae"\t"$4":"bs"-"be}'`
	r1=`echo $line | awk '{print $1}'`
	r2=`echo $line | awk '{print $2}'`
	echo $s faidx $r $r1
	$s faidx $r $r1 > $wd/$i/x1.fa
	echo $s faidx $r $r2
	$s faidx $r $r2 > $wd/$i/x2.fa
	cat $wd/$i/x1.fa $wd/$i/x2.fa > $wd/$i/x.fa

	line=`echo $rawline | awk -v f=1000 -v b=300 '{if ($1==0 || $1==2) {as=$3-f; ae=$3+b} else {as=$3-b; ae=$3+f}; if ($1==0 || $1==3) {bs=$5-b; be=$5+f} else {bs=$5-f; be=$5+b} ; if (as<=0) {as=1}; if (bs<=0) {bs=1}; print $2":"as"-"ae"\t"$4":"bs"-"be}'`
	r1=`echo $line | awk '{print $1}'`
	r2=`echo $line | awk '{print $2}'`

	
	echo $r | grep -i hg19 > /dev/null 
	if [ $? -eq 1 ]; then
		echo hg18
		#its hg18
		$s view -h $bam $r1 | python $g/scripts/filter.py $mean $stddev "" $wd/$i/out_r1.sam
		$s view $bam $r2 | python $g/scripts/filter.py $mean $stddev "" $wd/$i/out_r2.sam 
	else
		echo hg19 `echo $r1 | sed 's/chr//g'` `echo $r2 | sed 's/chr//g'`
		$s view -h $bam `echo $r1 | sed 's/chrx//g'` | python $g/scripts/filter.py $mean $stddev "" $wd/$i/out_r1.sam
		$s view $bam `echo $r2 | sed 's/chrx//g'` | python $g/scripts/filter.py $mean $stddev "" $wd/$i/out_r2.sam 
	fi
	cat $wd/$i/out_r1.sam $wd/$i/out_r2.sam | $s view -Sb -o $wd/$i/out.bam - 
	#$s view $bam $r1 | python $g/scripts/filter.py $mean $stddev "" $wd/out_r1.sam
	#$s view $bam $r2 | python $g/scripts/filter.py $mean $stddev "" $wd/out_r2.sam 
	#cat $wd/out_r1.sam $wd/out_r2.sam | awk '{print $1}' | sort | uniq > $wd/read_names
	#grep -f 
	#cat $wd/out_r1.sam $wd/out_r2.sam | $s view -Sb -o $wd/out.bam - 
	$s sort $wd/$i/out.bam $wd/$i/out_sorted
	$s view -h $wd/$i/out_sorted.bam > $wd/$i/out.sam
	java -jar $p/SamToFastq.jar I=$wd/$i/out.sam F=$wd/$i/f.fq F2=$wd/$i/r.fq
	python $g/scripts/trim_fastq.py $wd/$i/r.fq > $wd/$i/r_trim_.fq
	python $g/scripts/trim_fastq.py $wd/$i/f.fq > $wd/$i/f_trim_.fq
	python $g/scripts/remove_bad_quals.py $wd/$i/f_trim_.fq $wd/$i/r_trim_.fq $wd/$i/f_trim.fq $wd/$i/r_trim.fq
	cat $wd/$i/r_trim.fq | awk '{if (NR%4==1 || NR%4==2) {print $0}}' | sed 's/^@/>/g' > $wd/$i/r.fa
	cat $wd/$i/f_trim.fq | awk '{if (NR%4==1 || NR%4==2) {print $0}}' | sed 's/^@/>/g' > $wd/$i/f.fa
	#$v/velveth $wd/vel 51 -fasta -shortPaired -separate $wd/f.fa $wd/r.fa
	#$v/velvetg $wd/vel -min_pair_count 2 -ins_length $mean -scaffolding yes


	#assemble break point
	cat $wd/$i/f_trim.fq | awk '{if (NR%4==2) {print $0}}' > $wd/$i/fa
	cat $wd/$i/r_trim.fq | awk '{if (NR%4==2) {print $0}}' > $wd/$i/ra
	python $g/scripts/drop_short.py $wd/$i/fa $wd/$i/ra $wd/$i/f $wd/$i/r 30
	linesmall=`echo $rawline | awk -v f=1320 -v b=-275 '{if ($1==0 || $1==2) {as=$3-f; ae=$3+b} else {as=$3-b; ae=$3+f}; if ($1==0 || $1==3) {bs=$5-b; be=$5+f} else {bs=$5-f; be=$5+b} ; if (as<=0) {as=1}; if (bs<=0) {bs=1}; print $2":"as"-"ae"\t"$4":"bs"-"be}'`
	#linesmall=`echo $rawline | awk -v f=320 -v b=-275 '{if ($3=="-") {as=$2-b; ae=$2+f} else {as=$2-f; ae=$2+b}; if ($6=="-") {bs=$5-f; be=$5+b} else {bs=$5-b; be=$5+f} ; if (as<=0) {as=1}; if (bs<=0) {bs=1}; print $1":"as"-"ae"\t"$4":"bs"-"be}'`
	
	r1s=`echo $linesmall | awk '{print $1}'`
	r2s=`echo $linesmall | awk '{print $2}'`
	echo $r1s
	echo $r2s
	r1ss=`$s faidx $r $r1s | xargs echo | awk '{$1=""; print $0}' | sed 's/ //g'`
	rr1ss=`python $g/scripts/reverse.py $r1ss`
	r2ss=`$s faidx $r $r2s | xargs echo | awk '{$1=""; print $0}' | sed 's/ //g'`
	rr2ss=`python $g/scripts/reverse.py $r2ss`

#[misko@dupa05 minichr]$ cat /tmp/test | awk '{OFS="\t"; type=0; if ($3=="+") {if ($6=="+") {type=0} else {type=2} } else { if ($6=="+") {type=3} else {type=1} }; if ($1>$4 || ($1==$4 && $2>$5) ) { if (type<2) {type=1-type} ; fc=$4; fcc=$5; tc=$1; tcc=$2;} else {fc=$1; fcc=$2; tc=$4; tcc=$5;} ; print type,fc,fcc,tc,tcc}' | sort | uniq
#0       chr12   60846840        chr7    54413627
	typ=`echo $rawline | awk '{print $1}'`
	if [ $typ -eq 0 ] ; then
		$n $wd/$i/f $wd/$i/r $r1ss 25 13 $wd/$i/outa.fa 48 3000 > /dev/null 
		$n $wd/$i/f $wd/$i/r $rr2ss 25 13 $wd/$i/outb.fa 48 3000 > /dev/null
	elif [ $typ -eq 1 ]; then
		$n $wd/$i/f $wd/$i/r $rr1ss 25 13 $wd/$i/outa.fa 48 3000 > /dev/null
		$n $wd/$i/f $wd/$i/r $r2ss 25 13 $wd/$i/outb.fa 48 3000 > /dev/null
	elif [ $typ -eq 2 ]; then
		$n $wd/$i/f $wd/$i/r $r1ss 25 13 $wd/$i/outa.fa 48 3000 > /dev/null
		$n $wd/$i/f $wd/$i/r $r2ss 25 13 $wd/$i/outb.fa 48 3000 > /dev/null
	elif [ $typ -eq 3 ]; then
		$n $wd/$i/f $wd/$i/r $rr1ss 25 13 $wd/$i/outa.fa 48 3000 > /dev/null
		$n $wd/$i/f $wd/$i/r $rr2ss 25 13 $wd/$i/outb.fa 48 3000 > /dev/null
	fi

	$b -stepSize=5 -repMatch=2253 -minScore=0 -minIdentity=0 -maxIntron=0 -noHead -t=dna -q=dna $wd/$i/x.fa $wd/$i/outb.fa $wd/$i/outb.psl
	$b -stepSize=5 -repMatch=2253 -minScore=0 -minIdentity=0 -maxIntron=0 -noHead -t=dna -q=dna -out=axt $wd/$i/x.fa $wd/$i/outb.fa $wd/$i/outb.axt
	cat $wd/$i/outb.psl | python $g/scripts/resolve_bp.py $wd/$i/outb.fa $wd/$i/x1.fa $wd/$i/x2.fa > $wd/$i/results
	$b -stepSize=5 -repMatch=2253 -minScore=0 -minIdentity=0 -maxIntron=0 -noHead -t=dna -q=dna $wd/$i/x.fa $wd/$i/outa.fa $wd/$i/outa.psl
	$b -stepSize=5 -repMatch=2253 -minScore=0 -minIdentity=0 -maxIntron=0 -noHead -t=dna -q=dna -out=axt $wd/$i/x.fa $wd/$i/outa.fa $wd/$i/outa.axt
	cat $wd/$i/outa.psl | python $g/scripts/resolve_bp.py $wd/$i/outa.fa $wd/$i/x1.fa $wd/$i/x2.fa >> $wd/$i/results

	#rm $wd/$i/*.sam $wd/$i/*.bam
	#rm $wd/$i/*.fq $wd/$i/r $wd/$i/f $wd/$i/r.fa $wd/$i/f.fa
	

	i=`expr $i + 1`
done < $wd/uniq_links
