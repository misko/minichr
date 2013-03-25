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
#define MAX_STDDEV	23

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

class bp {
	public:
		unsigned int lefts;
		unsigned int rights;
		bp();
			
};

//some globals

double mean=0.0;
double stddev=0.0;

map < pos, bp > bps; 




pos::pos(int chr, unsigned int coord, bool strand) {
	this->chr=chr;
	this->coord=coord;
	this->strand=strand;
	this->marked=false;
}


bp::bp() {
	lefts=0;
	rights=0;
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


void update_bp(pos & a, pos & b) {

	if (a>b) {
		pos t = a;
		a=b;
		b=t;
	}

	if (a>b || a.chr!=b.chr) {
		cerr << " Failed pre condition " << endl;
		exit(1);
	}


	//find the nearest bp
	map<pos, bp>::iterator it= bps.lower_bound(a);
	while (it!=bps.end() && it->first.chr==a.chr && it->first.coord>=a.coord && it->first.coord<=b.coord) {
		if (a.marked) {
			it->second.lefts++;
		}
		if (b.marked) {
			it->second.rights++;
		}		
		it++;		
	}
		
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

			if (isize>(WEIRD_STDDEV*stddev+mean) || mate_chr!=my_chr) {
				//this is kinda weird
				return;
			}

			if (mate_pos < my_pos && mate_strand) {
				return; //cant figure out the rightmost part of mate
			}
			//cerr << my_chr << ":" << my_pos << " " << mate_chr << ":" << mate_pos << endl;

			pos my =  pos(my_chr,my_pos,my_strand);
			my.marked=true;
			pos mate = pos(mate_chr,mate_pos,mate_strand);


			#pragma omp critical 
			{
				update_bp(my,mate);
			}	
			
			
		}
	}
	return ;
}

void read_bps(char * bps_filename) {
	char const row_delim = '\n';
	char const field_delim = '\t';
	ifstream input( bps_filename );
	for (string row; getline(input, row, row_delim); ) {
		string schr ;
		int chr;
		unsigned int coord;
		istringstream ss(row);
		ss >> schr >> coord; 
		//cerr << schr << " " << chr << " " << coord << endl;
		chr=to_chr(schr.c_str());	
		pos p = pos(chr,coord,true);
		bps[p]=bp();
	}
	cerr << " Done reading bps ... " << endl;	
}

int main(int argc, char ** argv) {
	if (argc<4) {
		cerr << argv[0]  << " mean stddev bps" << endl;
		exit(1);
	}
	omp_set_num_threads(THREADS);
	mean=atof(argv[1]);
	stddev=atof(argv[2]);
	char * bps_filename=argv[3];


	if (mean<0.1 || stddev<0.1) {
		cerr << "invalid mean or stddev" << endl;
		exit(1);
	}

	cerr << "Using mean " << mean << " and stddev " << stddev << endl;

	read_bps(bps_filename);


	//got from online to read tab input


	unsigned int total_read=0;
	


	char const field_delim = '\t';
	vector< string > rows;
	rows.reserve(READ_SIZE);
	while (true) {
		//read in a million
		rows.clear();
		char buffer[1024];
		unsigned int read = 0;
		
		char row_buffer[1024*5];
		while(fgets(row_buffer,1024*5,stdin)!=NULL) {
			read++;
			rows.push_back(string(row_buffer));
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


	for (map<pos,bp>::iterator it=bps.begin(); it!=bps.end() ; it++) {
		char buffer[10]; 
		if (it->first.chr==23) {
			buffer[0]='X'; buffer[1]='\0';
		} else if (it->first.chr==24) {
			buffer[0]='Y'; buffer[1]='\0';
		} else {
			sprintf(buffer, "%d",it->first.chr);
		}
		cout << "chr" << buffer << "\t" << it->first.coord << "\t" << MAX(it->second.lefts,it->second.rights) << endl;		
	}	

	return 0;	

}
