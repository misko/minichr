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

#include "ssw_cpp.h"

using namespace std;
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define UNMAPPED        0x4
#define M_UNMAPPED      0x8
#define REVERSE         0x10
#define M_REVERSE       0x20

#define WEIRD_STDDEV    4
#define MAX_STDDEV      15

#define READ_SIZE	100000
#define NEAR	1000

class pos {
	public:
		bool marked;
		bool sharp;
		int chr;
		unsigned int coord;
		bool strand; // true is positive
		pos(int chr, unsigned int coord,bool strand);
		pos();
		bool operator<(const pos &other) const;
		bool operator>(const pos &other) const;
		bool operator==(const pos &other) const;
		bool operator!=(const pos &other) const;
		unsigned int operator-(const pos &other) const;
};

class cluster {
	public:
		pos b1;
		pos b2;
		cluster(pos,pos);
};




class alignment {
	public:
		pos inside;
		bool clip;
		alignment();
		alignment(pos inside, bool clip);
};

class cread {
	public:
		string name;
		int id; //into vector
		alignment a;
		int cid;
		cread(string name, int id, alignment a, int cid);

};


double mean, stddev;
char * ref[25];
unsigned int ref_sizes[25];

vector<cluster> clusters;
map<pos, set<  int > > cluster_pos;
map<string, vector< cread > > reads;

pos::pos(int chr, unsigned int coord, bool strand) {
	this->chr=chr;
	this->coord=coord;
	this->strand=strand;
	this->marked=false;
	this->sharp=false;
}


pos::pos() {
	chr=-1;
	coord=-1;
	strand=false;
	marked=false;
}

cluster::cluster(pos b1, pos b2) {
	this->b1 = b1;
	this->b2 = b2;
}
	
alignment::alignment() {
	this->inside=pos(0,0,'+');
}

alignment::alignment(pos inside, bool clip) {
	this->inside=inside;
	this->clip=clip;
}
		
cread::cread(string name , int id, alignment a, int cid) {
	this->name=name;
	this->id=id;
	this->a=a;
	this->cid=cid;
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

unsigned int pos::operator-(const pos &other) const {
	if (chr!=other.chr) {
		return -1;
	}
	if (coord>other.coord) {
		return coord-other.coord;
	}
	return other.coord-coord;
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


set<int> find_clusters_for_pos(pos & a ) {
	pos left_bound = pos(a.chr,MAX(stddev*MAX_STDDEV,a.coord)-stddev*MAX_STDDEV,true);
	pos right_bound = pos(a.chr,a.coord+stddev*MAX_STDDEV,true);

	map<pos, set<int> >::iterator left_it = cluster_pos.lower_bound(left_bound);	
	map<pos, set<int> >::iterator right_it = cluster_pos.upper_bound(right_bound);	

	set<int> found;

	//keep going!
	while (left_it!=right_it) {
		pos p = left_it->first;
		if (p.chr==a.chr && p.strand==a.strand) {
			if ( (a.coord+stddev*MAX_STDDEV > p.coord) || (p.coord+stddev*MAX_STDDEV > a.coord) ) {
				set<int> & s = left_it->second;
				for (set<int>::iterator sit=s.begin(); sit!=s.end(); sit++) {
					found.insert(*sit);	
				}
			}
		}
		left_it++;
	}

	return found;

}

int find_cluster(pos & a, pos & b) {

	if (a>b) {
		pos t = a;
		a=b;
		b=t;
	}

	if (a>b) {
		cerr << " Failed pre condition " << endl;
		exit(1);
	}

	set<int> foota = find_clusters_for_pos(a);
	set<int> footb = find_clusters_for_pos(b);
	
	//find how big the intersection is
	set<int> intersection;
	for (set<int>::iterator sit=foota.begin(); sit!=foota.end(); sit++) {
		for (set<int>::iterator ssit=footb.begin(); ssit!=footb.end(); ssit++) {
			if (*sit==*ssit) {
				intersection.insert(*sit);
			}
		}
	}
	
	if (intersection.size()>1) {
		cerr << "intersection bigger than 1!!! " << endl;
		for (set<int>::iterator sit=intersection.begin(); sit!=intersection.end(); sit++) {
			cerr << "\t" << *sit  << endl;
		}
		return -1;
	}

	if (intersection.size()==1) {
		return (*intersection.begin());
	}
	//return pair<pos,pos>(a,b);	
	return -1;

}



void read_ref(char * filename) {
	

	size_t length=0;
	char * buffer = (char*)malloc(1024*1024*1024);
	if (buffer==NULL) {
		cerr << "ERROR IN AMLLOC" << endl;
		exit(1);
	}


	FILE * fptr = fopen(filename,"r");
	if (fptr==NULL) {
		cerr << "Failed to read reference " << endl;
		exit(1);
	}
	char line_buffer[1024*5];
	string name="";
	while(fgets(line_buffer,1024*5,fptr)!=NULL) {
		line_buffer[strlen(line_buffer)-1]='\0';
		if (line_buffer[0]=='>') {
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
				//cerr << chr << endl;
			}
			//setup the new
			length=0;	
			name=string(line_buffer);
		} else {
			//add to current contig
			const char * s = line_buffer;
			unsigned int l = strlen(line_buffer);
			char buffer2[l+1];
			for (unsigned int i=0; i<l; i++) {
				buffer2[i]=tolower(s[i]);
			}
			buffer2[l]='\0'; //just in case make sure its legit
			for (unsigned int i=0; i<l; i++) {
				buffer[length+i]=toupper(buffer2[i]);
			}
			//memcpy(buffer+length,buffer2,l);
			length+=l;
		}
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
		//cerr << chr << endl;
	}



}


string reverse_comp(string s) {
	char buffer[s.length()+1];
	const char * c = s.c_str();
	unsigned int l = s.length();
	for (int i=0; i<l; i++) {
		switch (c[i]) {
			case 'A':
				buffer[l-1-i]='T';
				break;
			case 'T':
				buffer[l-1-i]='A';
				break;
			case 'G':
				buffer[l-1-i]='C';
				break;
			case 'C':
				buffer[l-1-i]='G';
				break;
			case 'N':
				buffer[l-1-i]='N';
				break;
			default:
				cerr << "FAILDE REV COMP " << endl;
				exit(1);
		}
	}
	buffer[l]='\0';
	return string(buffer);
}

int align(pos & align_near , string & squery) {
	cout << "ALIGN NEAR POS " << align_near.chr << ":" << align_near.coord << endl;
	pos start = pos(align_near.chr,MAX(NEAR,align_near.coord)-NEAR,'+');
	pos end = pos(align_near.chr,MIN(ref_sizes[align_near.chr-1],align_near.coord+NEAR),'+');


	if (end.coord<start.coord) {
		cerr << "Error " << endl;
		exit(1);
	}
	char cref[3000];
	memcpy(cref,ref[start.chr-1]+start.coord-1,end.coord-start.coord+1);
	cref[end.coord-start.coord+1]='\0';

	string sref=string(cref);

	StripedSmithWaterman::Aligner aligner(2,10,58,1);
	StripedSmithWaterman::Filter filter;
	StripedSmithWaterman::Alignment swalignment;
	//cerr << "Ref: " << sref << endl;
	//cerr << "Qeury: " << squery << endl;
	aligner.Align(squery.c_str(), sref.c_str(), sref.size(), filter, &swalignment);
cout << "===== SSW result =====" << endl;
//cout << "REF: " << squery << endl;
//cout << "Q: " << sref << endl;
cout << "Best Smith-Waterman score:\t" << swalignment.sw_score << endl
//<< "Next-best Smith-Waterman score:\t" << swalignment.sw_score_next_best << endl
//<< "Reference start:\t" << swalignment.ref_begin << endl
//<< "Reference end:\t" << swalignment.ref_end << endl
<< "Query start:\t" << swalignment.query_begin << endl
<< "Query end:\t" << swalignment.query_end << endl
//<< "Next-best reference end:\t" << swalignment.ref_end_next_best << endl
<< "Number of mismatches:\t" << swalignment.mismatches << endl
<< "Cigar: " << swalignment.cigar_string << endl;
cout << "======================" << endl;
	return swalignment.sw_score;

}

void process_read(vector<string> v_row) {
	//lets try to make an alignment and find its cluster
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


			int cid = find_cluster(my,mate);
			if (cid!=-1) {
				cerr << cid << endl;
				alignment a = alignment(my,my.sharp);
				int id = reads[qname].size();
				cread r = cread(qname,id,a,cid);
			}
			/*#pragma omp critical 
			{
				update_bp(my,mate);
			}*/	
			
			
		} else if ( ((flags & M_UNMAPPED) == 0) && ((flags & UNMAPPED) !=0)  ) {
			string qname = v_row[0];
			int my_chr = to_chr(v_row[2].c_str());

			int mate_chr = my_chr;
			if (v_row[6].c_str()[0]!='=') {
				mate_chr = to_chr(v_row[6].c_str());	
			}
			unsigned long mate_pos = atol(v_row[7].c_str());
			bool mate_strand = ((flags & M_REVERSE)==0);

			//lets find out what clusters we can go to 
			pos mate = pos(mate_chr,mate_pos,mate_strand);
			set<int> cids = find_clusters_for_pos(mate);
			for (set<int>::iterator sit=cids.begin(); sit!=cids.end(); sit++) {
				//need to find out which b we hit, or just map to both for all
				cluster & c = clusters[*sit];
			
				unsigned int d1 = mate-c.b1;
				unsigned int d2 = mate-c.b2;
				if (d1==-1 && d2==-1) {
					cerr << "both wrong " << endl;
					continue;
				}
				if (d1==d2) {
					cerr << "unclear " << endl;
					continue;
				}
				
				cerr << *sit << endl;
				int b1s=0,b1rs=0,b2s=0,b2rs=0;
				string squery = string(v_row[9]);	
				string rsquery = reverse_comp(squery);
				b1s = align(c.b1,squery);
				b2s = align(c.b2,squery);
				b1rs = align(c.b1,rsquery);
				b2rs = align(c.b2,rsquery);
				int max=MAX(MAX(b1s,b2s),MAX(b1rs,b2rs));
				if (max>80) {
					cout << qname << "support" << endl;
					cout << c.b1.chr << ":" << c.b1.coord << " " << c.b2.chr << ":" << c.b2.coord << endl;
					cout << mate-c.b1 << " " << mate-c.b2 << endl;
					
				}
			}
		}
	}
	return ;
	
}

int main( int argc, char ** argv) {
	if (argc!=5) {
		fprintf(stderr,"%s mean std cluster_file hg19.fa\n", argv[0]);
		exit(1);
	}

	mean = atof(argv[1]);
	stddev = atof(argv[2]);
	char * cluster_filename = argv[3];
	char * ref_filename = argv[4];
	read_ref(ref_filename);

	//read in the clusters
	ifstream infile(cluster_filename);
        string line;
        string name="";
        while (std::getline(infile, line)) {
		string s_chra;
		int chra = 0;
		unsigned int coorda=0;
		string stranda;
		
		string s_chrb;
		int chrb=0;
		unsigned int coordb=0;
		string strandb;
	
		int support=0;

		istringstream ss(line);
	
		ss >> s_chra >> coorda >> stranda >> s_chrb >> coordb >> strandb  >> support;  


		pos b1 = pos(to_chr(s_chra.c_str()),coorda,'+');
		pos b2 = pos(to_chr(s_chrb.c_str()),coordb,'+');
		cluster c = cluster(b1,b2);

		//check if already exists
		bool exists=false;
		set<int> b1_foot;
		for (set<int>::iterator sit = cluster_pos[b1].begin(); sit!=cluster_pos[b1].end(); sit++) {
			b1_foot.insert(*sit);
		}
		for (set<int>::iterator sit = cluster_pos[b2].begin(); sit!=cluster_pos[b2].end(); sit++) {
			if (b1_foot.find(*sit)!=b1_foot.end()) {
				//found in both, skip this cluster
				exists=true;
			}
		}

		if (!exists) {
			//want to add it
			//get a new cluster id
			int cid = clusters.size();
			clusters.push_back(c);
			cluster_pos[b1].insert(cid);
			cluster_pos[b2].insert(cid);
		}
	}

	//read in the reads
	unsigned long total_read=0;
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

		for (unsigned int i=0; i<read; i++) {
			//process it
			vector<string> v_row;
			istringstream ss(rows[i]);
			for (string field; getline(ss, field, field_delim); ) {
				v_row.push_back(field);
			}
			process_read(v_row);
		}
	
		if (read!=READ_SIZE) {
			break;
		}

	}
	cerr << "CLEAN OUT" << endl;
	
	return 0;
}
