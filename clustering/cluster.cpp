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

#define FAR_AWAY	50000
#define MIN_CLUSTER_SIZE	3

#define UNMAPPED	0x4
#define M_UNMAPPED	0x8
#define REVERSE		0x10
#define M_REVERSE	0x20

#define WEIRD_STDDEV	4
#define MAX_STDDEV	23

#define MIN_SHARP	10
#define SHARP_BP	15

#define THREADS	6

#define READ_SIZE	100000


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
map< pos, unsigned int > sharps;




pos median_sharp(pos & p, bool * sharp) {
	unsigned int total=0;
	for (int i=-MIN(p.coord,SHARP_BP); i<SHARP_BP; i++) {
		pos z = pos(p.chr, p.coord - i,p.strand);
		if (sharps.find(z)!=sharps.end()) {
			total+=sharps[z];
		}
	}
	if (total>10) {
		unsigned int median=total/2;
		total=0;
		for (int i=-MIN(p.coord,SHARP_BP); i<SHARP_BP; i++) {
			pos z = pos(p.chr, p.coord - i,p.strand);
			if (sharps.find(z)!=sharps.end()) {
				total+=sharps[z];
				if (total>=median) {
					*sharp=true;
					return z;
				}
			}
		}
	}
	*sharp=false;
	return p;
}

void insert_sharp(pos & p) {
	//sharps[find_sharp(p)].insert(p);
	//cerr << "SHARP " << p.chr << ":" << p.coord << endl;
	sharps[p]++;
}

pos::pos(int chr, unsigned int coord, bool strand) {
	this->chr=chr;
	this->coord=coord;
	this->strand=strand;
	this->marked=false;
	this->sharp=false;
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

pos set_median(multiset<pos> s) {
	multiset<pos>::iterator bit = s.begin();
	for (unsigned int i=0; i<s.size()/2; i++) {
		bit++;
	}
	return *bit;
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


unsigned int cigar_len(const char * s, bool * sharp) {
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
				case 'S':
					*sharp=true;
				case 'I':
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


bool much_smaller(const pos & a, const pos & b) {
	if (a.chr<b.chr) {
		return true;
	}
	if (a.chr>b.chr) {
		return false;	
	}
	if (b.coord > a.coord + FAR_AWAY) {
		return true;
	}
	return false;
}

void clean_clusters(pos max_pos) {
	//go over all the clusters and take out the ones that can't get more support and are very weak
	set< pair< pos , pos > > to_remove;
	for (map<pos , map<pos, cluster> >::iterator mmit=clusters.begin(); mmit!=clusters.end(); mmit++) {
		const pos & first_pos = mmit->first;
		bool first_smaller = much_smaller(first_pos,max_pos);
		for (map<pos, cluster>::iterator mit=clusters[first_pos].begin(); mit!=clusters[first_pos].end(); mit++) {
			const pos & second_pos = mit->first;
			bool second_smaller = much_smaller(second_pos,max_pos);
			cluster & c = clusters[first_pos][second_pos];
			if (first_smaller && second_smaller && c.lefts.size()<MIN_CLUSTER_SIZE && c.rights.size()<MIN_CLUSTER_SIZE) {
				//trim this cluster
				to_remove.insert(pair<pos,pos>(first_pos,second_pos));
			}
		}	
	}
	
	//cerr << "Preparing to remove " << to_remove.size() << " clusters" << endl;
	for (set< pair<pos, pos> >::iterator sit=to_remove.begin(); sit!=to_remove.end(); sit++) {
		clusters[sit->first].erase(sit->second);
		if (clusters[sit->first].size()==0) {
			clusters.erase(sit->first);
		}
	}
		
}

void update_cluster(pos a, pos b) {

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
	cluster c = cluster(a.strand,b.strand);
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


pos process_read(vector<string> & v_row) {
	//ok now we have a line lets check it out
	pos min_pos = pos(0,0,true);
	if (v_row[0].c_str()[0]=='@') {
		//its header
		//cerr << row << endl;
	} else {
		int flags = atoi(v_row[1].c_str());
		//check if both are mapped
		if ((flags & UNMAPPED)==0) {
			string qname = v_row[0];
			int my_chr = to_chr(v_row[2].c_str());
			unsigned long my_pos = atol(v_row[3].c_str());
			bool my_strand = ((flags & REVERSE)==0);

			pos my =  pos(my_chr,my_pos,my_strand);
			min_pos=my;

			my.marked=true;
			string my_cigar = v_row[5];
			unsigned int c_len=cigar_len(my_cigar.c_str(),&(my.sharp));

			if ((flags & (UNMAPPED + M_UNMAPPED))==0) {

				int mate_chr = my_chr;
				if (v_row[6].c_str()[0]!='=') {
					mate_chr = to_chr(v_row[6].c_str());	
				}
				unsigned long mate_pos = atol(v_row[7].c_str());
				bool mate_strand = ((flags & M_REVERSE)==0);
			
				double isize=1000*(stddev+mean+10); //TODO: HARD THRESHOLD
				if (mate_chr==my_chr) {
					isize=atof(v_row[8].c_str());
					if (isize<0) {
						isize=-isize;
					}
				}
				
				pos mate = pos(mate_chr,mate_pos,mate_strand);

				if (min_pos>mate) {
					min_pos=mate;
				}

				if (isize<(WEIRD_STDDEV*stddev+mean)) {
					//this is kinda normal
					if (my.sharp) {
						if (!my.strand) {
							my.coord+=c_len;
						}
						#pragma omp critical 
						{
							insert_sharp(my);
						}
					}
					
					return min_pos;
				}


				if (my.strand) {
					my.coord+=c_len;
					if (my.sharp) {
						#pragma omp critical 
						{
							insert_sharp(my);
						}

					}					
				}

				//cerr << my_chr << ":" << my_pos << " " << mate_chr << ":" << mate_pos << endl;

				//pos mate = pos(mate_chr,mate_pos,mate_strand);


				#pragma omp critical 
				{
					update_cluster(my,mate);
				}	
				
				
			} else {

			}
			
		}
	}
	return min_pos;
}

int main(int argc, char ** argv) {
	//assume input is sorted
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
	pos max_pos(0,0,true);
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

		/*for (string row; getline(cin, row, row_delim); ) {
			read++;
			rows.push_back(row);
			if (read==READ_SIZE) {
				break;
			}
		}*/

		total_read+=read;
		cerr << "\r" << "read: " << total_read << "     "; 

		//now process the reads
		#pragma omp parallel
		{
			unsigned int threads = omp_get_num_threads();
			unsigned int thread_id = omp_get_thread_num();
			pos t_max_pos=pos(0,0,true);
			for (unsigned int i=0; i<read; i++) {
				if (i%threads==thread_id) {
					//process it
					vector<string> v_row;
					istringstream ss(rows[i]);
					for (string field; getline(ss, field, field_delim); ) {
						v_row.push_back(field);
					}
					pos current_min = process_read(v_row);
					if (t_max_pos.chr==0 || current_min>t_max_pos) {
						t_max_pos=current_min;
					}
					
				}
			}
			#pragma omp critical
			{
				if (t_max_pos.chr!=0 && (max_pos.chr==0 || max_pos<t_max_pos)) {
					max_pos=t_max_pos;
				}
			}
		}

		if (max_pos.chr!=0) {
			clean_clusters(max_pos);
		}

		cerr << max_pos.chr << " : " << max_pos.coord << endl;
	
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
			if ( c.lefts.size()==0 || c.rights.size()==0) {
				continue;
			}

		
			pos left_bound;
			if (c.left_strand) {
				left_bound=set_max(c.lefts);
			} else {
				left_bound=set_min(c.lefts);
			}
			bool s_left=false;
			left_bound=median_sharp(left_bound,&s_left);

			pos right_bound;
			if (c.right_strand) {
				right_bound=set_max(c.rights);
			} else {
				right_bound=set_min(c.rights);
			}
			bool s_right=false;
			right_bound=median_sharp(right_bound,&s_right);



			//get the type
			int type=5;
			if (c.left_strand == !c.right_strand) {
				if (c.left_strand) {
					cout << "0\t" << left_bound.chr << ":" << left_bound.coord << "\t";
					cout << right_bound.chr << ":" << right_bound.coord << "\t" << c.lefts.size() << "\t";
					cout << (s_left ? "*" : "-" ) << "/" << (s_right ? "*" : "-" ) << endl;
	
					cout << "1\t" << right_bound.chr << ":" << right_bound.coord << "\t";
					cout << left_bound.chr << ":" << left_bound.coord << "\t" << c.lefts.size() << "\t";
					cout << (s_right ? "*" : "-" ) << "/" << (s_left ? "*" : "-" ) << endl;
				} else {
					cout << "1\t" << left_bound.chr << ":" << left_bound.coord << "\t";
					cout << right_bound.chr << ":" << right_bound.coord << "\t" << c.lefts.size() << "\t";
					cout << (s_left ? "*" : "-" ) << "/" << (s_right ? "*" : "-" ) << endl;
	
					cout << "0\t" << right_bound.chr << ":" << right_bound.coord << "\t";
					cout << left_bound.chr << ":" << left_bound.coord  << "\t" << c.lefts.size() << "\t";
					cout << (s_right ? "*" : "-" ) << "/" << (s_left ? "*" : "-" ) << endl;
				}
			} else {
				if (c.left_strand) {
					cout << "2\t" << left_bound.chr << ":" << left_bound.coord << "\t";
					cout << right_bound.chr << ":" << right_bound.coord << "\t" << c.lefts.size() << "\t";
					cout << (s_left ? "*" : "-" ) << "/" << (s_right ? "*" : "-" ) << endl;
	
					cout << "2\t" << right_bound.chr << ":" << right_bound.coord << "\t";
					cout << left_bound.chr << ":" << left_bound.coord << "\t" << c.lefts.size() << "\t";
					cout << (s_right ? "*" : "-" ) << "/" << (s_left ? "*" : "-" ) << endl;
				} else {
					cout << "3\t" << left_bound.chr << ":" << left_bound.coord << "\t";
					cout << right_bound.chr << ":" << right_bound.coord << "\t" << c.lefts.size() << "\t";
					cout << (s_left ? "*" : "-" ) << "/" << (s_right ? "*" : "-" ) << endl;

					cout << "3\t" << right_bound.chr << ":" << right_bound.coord << "\t";
					cout << left_bound.chr << ":" << left_bound.coord << "\t" << c.lefts.size() << "\t";
					cout << (s_right ? "*" : "-" ) << "/" << (s_left ? "*" : "-" ) << endl;
				}

			} 

			
		}
	}
	

	return 0;	

}
