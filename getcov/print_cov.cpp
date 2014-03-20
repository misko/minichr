#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
using namespace std;

unsigned short get_chr(const char * s) {
	if (s[0]=='C' || s[0]=='c') {
		s+=3;
	}
	char tmp[10];
	int i=0; 
	for (i=0; i<10 && isalnum(s[i]); i++) {
		tmp[i]=s[i];
	}
	tmp[i]='\0';
	cout << tmp << endl;	
	if (strlen(tmp)>0) {
		if (tmp[0]=='x' || tmp[0]=='X') {
			return 23;
		}
		if (tmp[0]=='y' || tmp[0]=='Y') {
			return 24;
		}
		if (tmp[0]=='m' || tmp[0]=='M') {
			return 25;
		}
		return atoi(tmp);
	} else {
		return 0;
	}
}

void read_cov(char * filename, bool normal) {
	cerr << "Reading coverage from file " << filename << endl;

	FILE * fptr = fopen(filename,"r");
	if (fptr==NULL) {
		fprintf(stderr, "Failed to open file %s\n",filename);
		exit(1);
	}
	
	//find the size of one entry
	size_t soe = sizeof(unsigned short)+sizeof(unsigned int)+sizeof(unsigned short);
	size_t chunk = 30737418240L;
	size_t size_so_far = 0;
	char * buffer = (char*) malloc(chunk);
	if (buffer==NULL) {
		cerr << " FALLED TO MALLOC " << endl;
		exit(1);
	}
	
	
	//the real read loop
	while (!feof(fptr)) {
		size_t read = fread(buffer+size_so_far,1,chunk,fptr);
		size_so_far+=read;
		cerr << "Read so far " << size_so_far << endl;

		if (read==chunk) {
			cerr << "REALLOC" << endl;
			exit(1);
			buffer=(char*)realloc(buffer,size_so_far+chunk);
			if (buffer==NULL) {
				cerr << " FALLED TO REALLOC " << endl;
				exit(1);
			}
		}
	}
	
	//lets get a buffer to fit the file	
	cout << "Done reading file  " << size_so_far <<  endl;

	size_t sz = size_so_far;	

	unsigned int entries = sz/soe;
	unsigned long total_coverage=0;

	cout << "Finding the average" << endl;

	unsigned short chr=0;
	unsigned int coord=0;
	unsigned short cov=0;
	
	for (unsigned int i=0; i<entries; i++) {
		char* base = buffer+i*soe;
		chr=*((unsigned short *)base);
		if (chr==0) {
			continue;
		}
		base+=sizeof(unsigned short);
		coord=*((unsigned int *)base);
		base+=sizeof(unsigned int);
		cov=*((unsigned short *)base);

		total_coverage+=cov;
	}


	double average=((double)total_coverage)/entries;

	double sum=0;

	cout << "Average is " << average << " , now finding stddev" << endl;	
	for (unsigned int i=0; i<entries; i++) {
		char* base = buffer+i*soe;
		chr=*((unsigned short *)base);
		if (chr==0) {
			continue;
		}
		base+=sizeof(unsigned short);
		coord=*((unsigned int *)base);
		base+=sizeof(unsigned int);
		cov=*((unsigned short *)base);


		sum+=(cov-average)*(cov-average)/entries;
	}


	double stddev=sqrt(sum);
	
	cout << "Standard deviation is " << stddev << endl;

	


	cerr << "total: " << total_coverage << endl;
	free(buffer);
}


char ** read_fasta(char * filename) {

	char ** fsta = (char**)malloc(sizeof(char*)*26);
	if (fsta==NULL) {
		cerr << "NOT GOOD" << endl;
		exit(1);
	}

	FILE * f = fopen(filename,"r");
	if (f==NULL) {
		cerr << " NOT GOOD 2 " << endl;
	}


	char * ref = (char*)malloc(sizeof(char)*20000000000L);
	if (ref==NULL) {
		fprintf(stderr, "ERROR\n");
		exit(1);
	}

	size_t ret = fread(ref,1,20000000000L,f);
	cout << "Read " << ret << " from ref" << endl;


	char * buffer = (char*)malloc(sizeof(char)*1000000000L);
	size_t i =0; 
	size_t j =0;
	buffer[j]='\0';
	int ichr=0;
	bool in_header=false;
	for (i=0; i<ret ; i++) {
		if (in_header || ref[i]=='\n') {
			if (in_header && ref[i]=='\n') {
				in_header=false;
			}
			continue;
		} else if (ref[i]=='>') {
			cerr << "HEADER " << get_chr(ref+i+1) << endl;
			ichr=get_chr(ref+i+1);
			in_header=true;
			if (ichr!=0) {
				fsta[ichr-1]=(char*)malloc(sizeof(char)*(j+1));
				if (fsta[ichr-1]==NULL) {
					cerr << "MAJOR ERROR" << endl;
					exit(1);
				}
				memcpy(fsta[ichr-1],buffer,(j+1)*sizeof(char));
				j=0;
				cout << "READ " << ichr << " " << i << endl;
			}
		} else {
			buffer[j++]=ref[i];
		}
	}
			if (ichr!=0) {
				fsta[ichr-1]=(char*)malloc(sizeof(char)*(i+1));
				if (fsta[ichr-1]==NULL) {
					cerr << "MAJOR ERROR" << endl;
					exit(1);
				}
				memcpy(fsta[ichr-1],buffer,(i+1)*sizeof(char));
				i=0;
				cout << "READ " << ichr << " " << i << endl;
			}
	
}

int main (int argc, char ** argv) {
	if (argc!=3) {
		cerr<<argv[0]<<" in_file ref"<<endl;
		return 0;
	}

	char * coverage_filename=argv[1];
	char * fasta_filename=argv[2];

	read_fasta(fasta_filename);
	read_cov(coverage_filename,false);

}

//chr1    10000   N       3       AA^:A   D@+
