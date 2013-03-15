#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string.h>
#include <string>
#include <map>
#include <set>

using namespace std;
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define UNMAPPED	0x4
#define M_UNMAPPED	0x8
#define REVERSE		0x10
#define M_REVERSE	0x20

#define WEIRD_STDDEV	4
#define MAX_STDDEV	15

#define SIZE_RES	1
#define GC_RES	1

class pos {
	public:
		int chr;
		unsigned int coord;
		bool strand; // true is positive
		pos(int chr, unsigned int coord,bool strand);
		pos();
		bool operator<(const pos &other) const;
		bool operator>(const pos &other) const;
		bool operator==(const pos &other) const;
		bool operator!=(const pos &other) const;
};


//some globals

double mean=0.0;
double stddev=0.0;
unsigned int total=0;
map< int, unsigned int > size_totals;
map< string , pos > mappings;
map< int , map< int , unsigned int > > gcs;
char * ref[25];
size_t ref_sizes[25];
int max_size=0;

pos::pos(int chr, unsigned int coord, bool strand) {
	this->chr=chr;
	this->coord=coord;
	this->strand=strand;
}



pos::pos() {
	chr=-1;
	coord=-1;
	strand=false;
}

bool pos::operator<(const pos &other) const {
	if (chr<other.chr) {
		return true;
	} else if (other.chr==chr) {
		return coord<other.coord;
	} else {
		return false;
	}
}

bool pos::operator>(const pos &other) const {
	if (chr>other.chr) {
		return true;
	} else if (other.chr==chr) {
		return coord>other.coord;
	} else {
		return false;
	}
}

bool pos::operator==(const pos &other) const {
	if (chr!=other.chr || coord!=other.coord || strand!=other.strand) {
		return false;
	} 
	return true;

}


bool pos::operator!=(const pos &other) const {
	if (*this==other) {
		return false;
	}
	return true;

}



pos set_max(multiset<pos> s) {
	if (s.size()==0) {
		cerr << " INVALID SET " << endl;
		exit(1);
	}
	multiset<pos>::iterator it = s.begin();
	pos x = *it;
	it++;
	while (it!=s.end()) {
		if ((*it)>x) {
			x=*it;
		}
		it++;
	}
	return x;
}

pos set_min(multiset<pos> s) {
	if (s.size()==0) {
		cerr << " INVALID SET " << endl;
		exit(1);
	}
	multiset<pos>::iterator it = s.begin();
	pos x = *it;
	it++;
	while (it!=s.end()) {
		if ((*it)<x) {
			x=*it;
		}
		it++;
	}
	return x;
}

//take a chr name and give back the int
int to_chr(const char * s) {
	char buff[1024]="";
	unsigned int i=0;
	for (; i<strlen(s); i++) {
		buff[i]=tolower(s[i]);
	}
	buff[i]='\0';
	char * p = buff;
	if (i>3 && buff[0]=='c' && buff[1]=='h' && buff[2]=='r') {
		p=buff+3;
	}
	if (p[0]=='x') {
		return 23;
	}
	if (p[0]=='y') {
		return 24;
	}
	if (p[0]=='m') {
		return 25;
	}
	return atoi(p);
}



void read_ref(char * filename) {
	

	size_t length=0;
	char * buffer = (char*)malloc(1024*1024*1024);
	if (buffer==NULL) {
		cerr << "ERROR IN AMLLOC" << endl;
		exit(1);
	}

	ifstream infile(filename);
	string line;	
	string name="";
	while (std::getline(infile, line)) {
		if (line.c_str()[0]=='>') {
			//dump the last contig
			if (name.length()>0) {
				int chr= to_chr(name.c_str()+1);
				if (chr==0) {
					cerr << "Failed to read fasta" << endl;
					exit(1);
				}
				char * ref_chr = (char*)malloc(length);
				if (ref_chr==NULL) {
					cerr << "Failed to amlloc" << endl;
					exit(1);
				}
				memcpy(ref_chr,buffer,length);
				ref[chr-1]=ref_chr;
				ref_sizes[chr-1]=length;
				cerr << chr << endl;
			}
			//setup the new
			length=0;	
			name=line;
		} else {
			//add to current contig
			const char * s = line.c_str();
			char buffer2[line.length()+1];
			for (unsigned int i=0; i<line.length(); i++) {
				buffer2[i]=tolower(s[i]);
			}
			buffer2[line.length()]='\0'; //just in case make sure its legit
			memcpy(buffer+length,buffer2,line.length());
			length+=line.length();
		}
		//cerr << "X" << line << "Y" << endl;
	}
	if (name.length()>0) {
		int chr= to_chr(name.c_str()+1);
		if (chr==0) {
			cerr << "Failed to read fasta" << endl;
			exit(1);
		}
		char * ref_chr = (char*)malloc(length);
		if (ref_chr==NULL) {
			cerr << "Failed to amlloc" << endl;
			exit(1);
		}
		memcpy(ref_chr,buffer,length);
		ref[chr-1]=ref_chr;
		ref_sizes[chr-1]=length;
		cerr << chr << endl;
	}

}


int get_gc(pos left, pos right) {
	if (left.chr!=right.chr) {
		cerr << " wrong chr! "  << endl;
		exit(1);
	}
	int n=0;
	int gc=0;
	int at=0;
	for (unsigned int i=MAX(1,left.coord)-1; i<MIN(ref_sizes[left.chr-1],right.coord); i++) {
		char c = ref[left.chr-1][i];
		switch (c) {
			case 'a':
			case 't':
				at++;
				break;
			case 'c':
			case 'g':
				gc++;
				break;
			case 'n':
				n++;
				break;
			default:
				cerr << "NO CASE FOR " << c << endl;
				exit(1);
		}
	} 	
	return gc;
}


int main(int argc, char ** argv) {
	if (argc<4) {
		cerr << argv[0]  << " hg19.fa mean stddev" << endl;
		exit(1);
	}

	char * hg_filename=argv[1];
	mean=atof(argv[2]);
	stddev=atof(argv[3]);

	if (mean<0.1 || stddev<0.1) {
		cerr << "invalid mean or stddev" << endl;
		exit(1);
	}

	cerr << "Using mean " << mean << " and stddev " << stddev << endl;

	read_ref(hg_filename);
	
	//got from online to read tab input
	char const row_delim = '\n';
	char const field_delim = '\t';
	for (string row; getline(cin, row, row_delim); ) {
		vector<string> v_row;
		istringstream ss(row);
		for (string field; getline(ss, field, field_delim); ) {
			v_row.push_back(field);
		}
		//ok now we have a line lets check it out
		if (v_row[0].c_str()[0]=='@') {
			//its header
			//cerr << row << endl;
		} else {
			int flags = atoi(v_row[1].c_str());
			//check if both are mapped
			if ((flags & (UNMAPPED + M_UNMAPPED))==0) {
				string qname = v_row[0];
				int my_chr = to_chr(v_row[2].c_str());
				unsigned long my_pos = atol(v_row[3].c_str());
				bool my_strand = ((flags & REVERSE)==0);

				int mate_chr = my_chr;
				if (v_row[6].c_str()[0]!='=') {
					mate_chr = to_chr(v_row[6].c_str());	
				}
				unsigned long mate_pos = atol(v_row[7].c_str());
				bool mate_strand = ((flags & M_REVERSE)==0);
			
				double isize=1000*mean; //TODO: HARD THRESHOLD
				if (mate_chr==my_chr) {
					isize=atof(v_row[8].c_str());
				}

				if (isize>=(WEIRD_STDDEV*stddev+mean)) {
					//this is kinda weird
					continue;
				}
	

				//cerr << my_chr << ":" << my_pos << (my_strand ?  "+" : "-" )  << " " << mate_chr << ":" << mate_pos  << (mate_strand ? "+" : "-" ) << endl;

				if (mappings.find(qname)!=mappings.end()) {
					//cerr << "FOUND " << endl;
					//lets add it to the clusters
					pos my = pos(my_chr,my_pos,my_strand);
					pos mate = mappings[qname];
					mappings.erase(qname);	
					
					//process the read
					if (isize<0) {
						isize=-isize;
					}	
	
					pos left = my;
					pos right = mate;
					if (left>right) {
						left=mate;
						right=my;
					}
						
					if (left.strand && !right.strand) {
						//the only case we count it
						/*
						if (isize>max_size) {
							max_size=isize;
						}
						int gc = get_gc(left,pos(right.chr,left.coord+isize,right.strand));
						gcs[SIZE_RES*((int)(isize/SIZE_RES))][GC_RES*((int)(gc/GC_RES))]++;
						total++;
						size_totals[SIZE_RES*((int)isize/SIZE_RES)]++; */
						max_size=400;
						unsigned int center = left.coord + isize/2;
						pos xleft = pos(left.chr, MAX(200,center)-200, left.strand); 
						pos xright = pos(right.chr, center+200, right.strand); 
						int gc = get_gc(xleft,xright);
						gcs[400][gc]++;
						total++;
						size_totals[400]++; 
						
						//cerr << "HAS GC " << gc << endl;
					}
				} else {
					//need to add it to mappings
					//cerr << "ADDING " << qname << endl;
					mappings[qname]=pos(my_chr,my_pos,my_strand);
				}
				
			}
		}
	}

	

	//print out gc matrix
	/*cout << "SIZE\t";
	for (int j=0; j*GC_RES<=max_size; j++) {
		cout << j*GC_RES << "\t" ;
	}
	cout << endl;

	for (int i=1; i*SIZE_RES<=max_size; i++) {
		if ((total/20)<size_totals[i*SIZE_RES]) {
			cout << i*SIZE_RES << "\t";
			for (int j=0; j*GC_RES<=max_size; j++) {
				cout << gcs[i*SIZE_RES][j*GC_RES] << "\t" ;
			}
			cout << endl;
		}
	}

	cout << endl << endl;*/


	//print out fractions	
	for (int i=1; i*SIZE_RES<=max_size; i++) {
		if ((total/20)<size_totals[i*SIZE_RES]) {
			cout << i*SIZE_RES << "\t";
			for (int j=0; j*GC_RES<=i*SIZE_RES; j++) {
				cout << ((double)j*GC_RES)/i*SIZE_RES << "\t" ;
			}
			cout << endl;

			cout << i*SIZE_RES << "\t";
			for (int j=0; j*GC_RES<=i*SIZE_RES; j++) {
				cout << gcs[i*SIZE_RES][j*GC_RES] << "\t" ;
			}
			cout << endl;
		}
	}

	cout << endl << endl;

	return 0;	

}
