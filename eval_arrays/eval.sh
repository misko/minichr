fn=$1
len=$2
n=$3
sp=$4

/bin/hostname
rm /dupa-filer/misko/tcga_cn_array/fits/`basename $1`_l${len}_n${n}_sp${sp}
rm /dupa-filer/misko/tcga_cn_array/fits/`basename $1`_l${len}_n${n}_sp${sp}.gz

while read line; do 
	/usr/bin/pypy /dupa-filer/misko/tcga_cn_array/cfy.py $line $len $n $sp	 >> /dupa-filer/misko/tcga_cn_array/fits/`basename $1`_l${len}_n${n}_sp${sp}
done < $1

/bin/gzip /dupa-filer/misko/tcga_cn_array/fits/`basename $1`_l${len}_n${n}_sp${sp}

/bin/hostname
