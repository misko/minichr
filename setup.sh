killall pigz; rm -f inT inN; mkfifo inT inN;

/filer/misko/mini_chr/git/minichr/pigz-2.2.5/unpigz -c -p 12 /dev/shm/DIPG29T.cov.gz > inT &
/filer/misko/mini_chr/git/minichr/pigz-2.2.5/unpigz -c -p 12 /dev/shm/DIPG29N_replacement.cov.gz > inN & 
