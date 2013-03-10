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
#define REVERSE		0x1
#define M_REVERSE	0x2

#define WEIRD_STDDEV	4
#define MAX_STDDEV	15


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

map< string , pos > mappings;
map< pos, map < pos, cluster > > clusters; 

pos::pos(int chr, unsigned int coord, bool strand) {
	this->chr=chr;
	this->coord=coord;
	this->strand=strand;
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

pair<pos, pos> find_cluster(pos & a, pos & b) {

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
								return pair<pos,pos>(p,xp);
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
	clusters[a][b]=cluster(a.strand,b.strand);
	
	return pair<pos,pos>(a,b);	


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




int main(int argc, char ** argv) {
	if (argc<3) {
		cerr << argv[0]  << " mean stddev" << endl;
		exit(1);
	}

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
	for (string row; getline(cin, row, row_delim); ) {
		vector<string> v_row;
		istringstream ss(row);
		for (string field; getline(ss, field, field_delim); ) {
			v_row.push_back(field);
		}
		//ok now we have a line lets check it out
		if (v_row[0].c_str()[0]=='@') {
			//its header
			cerr << row << endl;
		} else {
			int flags = atoi(v_row[1].c_str());
			//check if both are mapped
			if ((flags & (UNMAPPED + M_UNMAPPED))==0) {
				string qname = v_row[0];
				int my_chr = to_chr(v_row[2].c_str());
				unsigned long my_pos = atol(v_row[3].c_str());
				bool my_strand = !(flags & REVERSE);

				int mate_chr = my_chr;
				if (v_row[6].c_str()[0]!='=') {
					mate_chr = to_chr(v_row[6].c_str());	
				}
				unsigned long mate_pos = atol(v_row[7].c_str());
				bool mate_strand = !(flags & M_REVERSE);
			
				double isize=1000*mean; //TODO: HARD THRESHOLD
				if (mate_chr==my_chr) {
					isize=atof(v_row[8].c_str());
				}

				if (isize<(WEIRD_STDDEV*stddev+mean)) {
					//this is kinda normal
					continue;
				}
	

				cerr << my_chr << ":" << my_pos << " " << mate_chr << ":" << mate_pos << endl;

				if (mappings.find(qname)!=mappings.end()) {
					cerr << "FOUND " << endl;
					//lets add it to the clusters
					pos my = pos(my_chr,my_pos,my_strand);
					pos mate = mappings[qname];
					mappings.erase(qname);	
					
					pair<pos,pos> key = find_cluster(my,mate);
					pos left = my;
					pos right = mate;
					if (my>mate) {
						left=mate;
						right=my;
					}

					clusters[key.first][key.second].lefts.insert(left);					
					clusters[key.first][key.second].rights.insert(right);					

				} else {
					//need to add it to mappings
					cerr << "ADDING " << qname << endl;
					mappings[qname]=pos(my_chr,my_pos,my_strand);
				}
				
			}
		}
	}

					cerr << "CLUSTER FOUND " << endl;
	//print out the clusters
	for (map<pos , map< pos, cluster> >::iterator it=clusters.begin(); it!=clusters.end(); it++) {
		pos posa = it->first;
		for (map<pos,cluster>::iterator xit=clusters[it->first].begin(); xit!=clusters[it->first].end(); xit++) {
			pos posb = xit->first;
			cluster c = clusters[posa][posb];
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
			cerr << "CLUSTER" << endl;	
			cerr << left_bound.chr << ":" << left_bound.coord << left_bound.strand << "\t" << right_bound.chr << ":" << right_bound.coord << right_bound.strand << " " << c.lefts.size() << endl;
		}
	}
	

	return 0;	

}
