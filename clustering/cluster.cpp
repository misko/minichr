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
#include <omp.h>

using namespace std;
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define UNMAPPED	0x4
#define M_UNMAPPED	0x8
#define REVERSE		0x10
#define M_REVERSE	0x20

#define WEIRD_STDDEV	4
#define MAX_STDDEV	15

#define THREADS	6

#define READ_SIZE	1000000

class pos {
	public:
		bool marked;
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

class cluster {
	public:
		multiset<pos> lefts;
		bool left_strand;
		multiset<pos> rights;
		bool right_strand;
		cluster(bool left_strand, bool right_strand);
		cluster();
			
};

//some globals

double mean=0.0;
double stddev=0.0;

map< pos, map < pos, cluster > > clusters; 

pos::pos(int chr, unsigned int coord, bool strand) {
	this->chr=chr;
	this->coord=coord;
	this->strand=strand;
	this->marked=false;
}


cluster::cluster(bool left_strand, bool right_strand) {
	this->left_strand=left_strand;
	this->right_strand=right_strand;
}

cluster::cluster() {

}

pos::pos() {
	chr=-1;
	coord=-1;
	strand=false;
	marked=false;
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


unsigned int cigar_len(const char * s) {
	unsigned int len=0;
	unsigned int xlen=0;
	for (int i=0; i<strlen(s); i++) {
		if (isdigit(s[i])) {
			//keep going
			xlen=xlen*10+(s[i]-48);
		} else {
			//ok lets process the op
			switch(s[i]) {
				case 'M':
				case 'D':
					len+=xlen;
					break;
				case 'I':
				case 'S':
					break;
				default:
					cerr << "Failed to handle cigar op " << s[i] << endl;
					exit(1);		
			}
			xlen=0;
		}
	}
	//cerr << "CIGAR " << s << " " << len << endl;
	return len;
}


void update_cluster(pos & a, pos & b) {

	if (a>b) {
		pos t = a;
		a=b;
		b=t;
	}

	if (a>b) {
		cerr << " Failed pre condition " << endl;
		exit(1);
	}

	pos left_bound = pos(a.chr,MAX(stddev*MAX_STDDEV,a.coord)-stddev*MAX_STDDEV,true);
	pos right_bound = pos(a.chr,a.coord+stddev*MAX_STDDEV,true);
	pos xleft_bound = pos(b.chr,MAX(stddev*MAX_STDDEV,b.coord)-stddev*MAX_STDDEV,true);
	pos xright_bound = pos(b.chr,b.coord+stddev*MAX_STDDEV,true);
	map<pos , map<pos, cluster> >::iterator left_it = clusters.lower_bound(left_bound);
	map<pos , map<pos, cluster> >::iterator right_it = clusters.upper_bound(right_bound);
	if(left_it!=right_it) {
		//keep going!
		while (left_it!=right_it) {
			pos p = left_it->first;
			if (p.chr==a.chr && p.strand==a.strand) {
				if ( (a.coord+stddev*MAX_STDDEV > p.coord) || (p.coord+stddev*MAX_STDDEV > a.coord) ) {
					//print have a match for this round
					map<pos, cluster>::iterator xleft_it = clusters[p].lower_bound(xleft_bound);
					map<pos, cluster>::iterator xright_it = clusters[p].upper_bound(xright_bound);
					while (xleft_it!=xright_it) {
						pos xp = xleft_it->first;
						if (xp.chr==b.chr && xp.strand==b.strand) {
							if ((b.coord+stddev*MAX_STDDEV > xp.coord) || (xp.coord+stddev*MAX_STDDEV > b.coord)) {
								//found a match!
								//return pair<pos,pos>(p,xp);
								if (a.marked) {
									(xleft_it->second).lefts.insert(a);
								}
								if (b.marked) {
									(xleft_it->second).rights.insert(b);
								}
								return;
							}
						}	
						xleft_it++;
					}	
				}	
			}
			left_it++;
		}
	}
	
	//make a new entry
	cluster c  =cluster(a.strand,b.strand);
	if (a.marked) {
		c.lefts.insert(a);
	}
	if (b.marked) {
		c.rights.insert(b);
	}
	clusters[a][b]=c;
	
	//return pair<pos,pos>(a,b);	
	return ;

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



#define P	7919 

unsigned int hash_string(const char * s) {
	unsigned int z=0;
	for (int i=0; i<strlen(s); i++) {
		z+=s[i]*P+z*(P+1);
	}
	//cerr << "HASH " << s << " -> " << z << endl;
	return z;
}


void process_read(vector<string> & v_row) {
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

			if (my_strand) {
				string my_cigar = v_row[5];
				my_pos+=cigar_len(my_cigar.c_str());
			}

			int mate_chr = my_chr;
			if (v_row[6].c_str()[0]!='=') {
				mate_chr = to_chr(v_row[6].c_str());	
			}
			unsigned long mate_pos = atol(v_row[7].c_str());
			bool mate_strand = ((flags & M_REVERSE)==0);
		
			double isize=1000*mean; //TODO: HARD THRESHOLD
			if (mate_chr==my_chr) {
				isize=atof(v_row[8].c_str());
				if (isize<0) {
					isize=-isize;
				}
			}

			if (isize<(WEIRD_STDDEV*stddev+mean)) {
				//this is kinda normal
				return;
			}


			//cerr << my_chr << ":" << my_pos << " " << mate_chr << ":" << mate_pos << endl;

			pos my =  pos(my_chr,my_pos,my_strand);
			my.marked=true;
			pos mate = pos(mate_chr,mate_pos,mate_strand);


			#pragma omp critical 
			{
				update_cluster(my,mate);
			}	
			
			
		}
	}
	return ;
}

int main(int argc, char ** argv) {
	if (argc<3) {
		cerr << argv[0]  << " mean stddev" << endl;
		exit(1);
	}
	omp_set_num_threads(THREADS);
	mean=atof(argv[1]);
	stddev=atof(argv[2]);

	if (mean<0.1 || stddev<0.1) {
		cerr << "invalid mean or stddev" << endl;
		exit(1);
	}

	cerr << "Using mean " << mean << " and stddev " << stddev << endl;


	//got from online to read tab input
	char const row_delim = '\n';
	char const field_delim = '\t';


	unsigned int total_read=0;
	


	vector< string > rows;
	rows.reserve(READ_SIZE);
	while (true) {
		//read in a million
		char buffer[1024];
		unsigned int read = 0;
		for (string row; getline(cin, row, row_delim); ) {
			read++;
			rows.push_back(row);
			if (read==READ_SIZE) {
				break;
			}
		}

		total_read+=read;
		cerr << "\r" << "read: " << total_read << "     "; 

		//now process the reads
		#pragma omp parallel
		{
			unsigned int threads = omp_get_num_threads();
			unsigned int thread_id = omp_get_thread_num();
			for (unsigned int i=0; i<read; i++) {
				if (i%threads==thread_id) {
					//process it
					vector<string> v_row;
					istringstream ss(rows[i]);
					for (string field; getline(ss, field, field_delim); ) {
						v_row.push_back(field);
					}
					process_read(v_row);
				}
			}
		}
	
		if (read!=READ_SIZE) {
			break;
		}

	}
	cerr << "CLEAN OUT" << endl;
	//print out the clusters
	for (map<pos , map< pos, cluster> >::iterator it=clusters.begin(); it!=clusters.end(); it++) {
		pos posa = it->first;
		for (map<pos,cluster>::iterator xit=clusters[it->first].begin(); xit!=clusters[it->first].end(); xit++) {
			pos posb = xit->first;
			cluster c = clusters[posa][posb];
			if (c.lefts.size()==0 || c.rights.size()==0) {
				continue;
			}
			pos left_bound;
			if (c.left_strand) {
				left_bound=set_max(c.lefts);
			} else {
				left_bound=set_min(c.lefts);
			}
			pos right_bound;
			if (c.right_strand) {
				right_bound=set_max(c.rights);
			} else {
				right_bound=set_min(c.rights);
			}

			//get the type
			int type=5;
			if (c.left_strand == !c.right_strand) {
				if (c.left_strand) {
					cout << "0\t" << left_bound.chr << ":" << left_bound.coord << left_bound.strand << "\t";
					cout << right_bound.chr << ":" << right_bound.coord << right_bound.strand << " " << c.lefts.size() << endl;
	
					cout << "1\t" << right_bound.chr << ":" << right_bound.coord << right_bound.strand << "\t";
					cout << left_bound.chr << ":" << left_bound.coord << left_bound.strand << " " << c.lefts.size() << endl;
				} else {
					cout << "1\t" << left_bound.chr << ":" << left_bound.coord << left_bound.strand << "\t";
					cout << right_bound.chr << ":" << right_bound.coord << right_bound.strand << " " << c.lefts.size() << endl;
	
					cout << "0\t" << right_bound.chr << ":" << right_bound.coord << right_bound.strand << "\t";
					cout << left_bound.chr << ":" << left_bound.coord << left_bound.strand << " " << c.lefts.size() << endl;
				}
			} else {
				if (c.left_strand) {
					cout << "2\t" << left_bound.chr << ":" << left_bound.coord << left_bound.strand << "\t";
					cout << right_bound.chr << ":" << right_bound.coord << right_bound.strand << " " << c.lefts.size() << endl;
	
					cout << "2\t" << right_bound.chr << ":" << right_bound.coord << right_bound.strand << "\t";
					cout << left_bound.chr << ":" << left_bound.coord << left_bound.strand << " " << c.lefts.size() << endl;
				} else {
					cout << "3\t" << left_bound.chr << ":" << left_bound.coord << left_bound.strand << "\t";
					cout << right_bound.chr << ":" << right_bound.coord << right_bound.strand << " " << c.lefts.size() << endl;

					cout << "3\t" << right_bound.chr << ":" << right_bound.coord << right_bound.strand << "\t";
					cout << left_bound.chr << ":" << left_bound.coord << left_bound.strand << " " << c.lefts.size() << endl;
				}

			} 
		}
	}
	

	return 0;	

}
