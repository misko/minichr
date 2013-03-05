#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <vector>
#include <set>
#include <map>
#include <fstream>
#include <sstream>
#include <iostream>
#include <zlib.h>
#include <omp.h>
#include <queue>
#include <math.h>

#define THREADS	32

#define MIN_FLOW 4
#define MAX_FLOW 30
#define MAX_FREE 20

#define ZERO 1e-5

#define MAX_EDGE_SIZE 2400
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

using namespace std;


double abs(double a) {
	if (a<0) {
		return -a;
	}
	return a;
}

class pos {
	public:
		int chr;
		unsigned int coord;
		pos(int chr, unsigned int coord);
		pos();
		bool operator<(const pos &other) const;
		bool operator>(const pos &other) const;
		bool operator==(const pos &other) const;
		bool operator!=(const pos &other) const;
};

class edge {
	public:
		pos posa;
		pos posb;
		edge(pos posa, pos posb);
		edge reverse();
		unsigned int length();
		bool operator<(const edge &other) const;
		bool operator>(const edge &other) const;
		bool operator==(const edge &other) const;
		bool is_forward();
};


class edge_info { 
	public:
		edge_info();
		unsigned int bp;
		double normal_coverage;
		double cancer_coverage;
		unsigned int type;
		unsigned int supporting;
};

class state {
	public:
		vector<double> scores;
		//store the path as vector and multiset
		vector<edge> gpath_vector;
		multiset<edge> gpath_set;
		//store the free edges used
		multiset<edge> fpath_set;
		bool bp_check;
		double score;
		int cp;
		double ncov,ccov,pcov;
		int max_flow;	

		//constructors
		state(edge e);
		state();
		state(state s, edge e);
		state(vector<edge> ve);
	
		//constructor helpers
		void prepend_edge(edge e);
		void append_edge(edge e);
		void add_edge(edge e, bool append);
		void add_edge_to_score(edge e);

		//successor function
		vector<state> children();

		//scoring
		void best_score();
		double score_with_flow(int flow);
		bool go_on();
		double tail_check();
		bool is_dup();
		unsigned int length();

		//opeartors
		bool operator<(const state &other) const;
		bool operator>(const state &other) const;

		//display
		string str();
};




class state_hash {
	public:
		unsigned int last_coord;
		unsigned int first_coord; 	
		double score;
		int cp;
		double ncov,ccov,pcov;
		state_hash(state s) {
			last_coord = first_coord =0;
			if (s.gpath_vector.size()>0) {
				last_coord = s.gpath_vector.back().posb.coord;
				first_coord = s.gpath_vector.front().posa.coord;
			}
			
			score=s.score;
			cp=s.cp;
			ncov=s.ncov;
			ccov=s.ccov;
			pcov=s.pcov;
		}
		//opeartors
		bool operator<(const state_hash &other) const {
			if (last_coord+first_coord<other.last_coord+other.first_coord) {
				return true;
			}
			return false;
		}
		bool operator==(const state_hash &other) const {
			return memcmp(this,&other,sizeof(state_hash))==0;
		}
		size_t operator()(const state_hash &other) const {
			return last_coord+first_coord;
		}
};


//some global variables
set<pos> bps;
map<edge, edge_info > edges;
map<pos, set<pos> > free_edges;
map<pos, double> cancer_pair_coverage;
map<pos, double> normal_pair_coverage;
unsigned int total_normal_pair_arcs;
unsigned int total_cancer_pair_arcs;
map<edge, state> edge_bests;


unsigned int bp_range=100000;

edge fake_edge=edge(pos(0,0),pos(0,0));

//
// The pos class
//


pos::pos(int chr, unsigned int coord) {
	this->chr=chr;
	this->coord=coord;
}

pos::pos() {
	chr=-1;
	coord=-1;
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
	if (chr!=other.chr || coord!=other.coord) {
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



//
// The edge class
//
//if supporting > 0 then its a free edge!
bool edge::operator<(const edge &other) const {
	if (posa<other.posa) {
		return true;
	} else if (posa>other.posa) {
		return false;
	} else {
		return posb<other.posb;
	}
}

bool edge::operator>(const edge &other) const {
	if (posa>other.posa) {
		return true;
	} else if (posa<other.posa) {
		return false;
	} else {
		return posb>other.posb;
	}
}

edge edge::reverse() {
	return edge(posb,posa);
}

unsigned int edge::length() {
	if (posa.chr==posb.chr) {
		if (posa.coord>posb.coord) {
			return posa.coord-posb.coord;
		} else {
			return posb.coord-posa.coord;
		}
	}
	return 0;
}

bool edge::operator==(const edge &other) const {
	return (posa==other.posa && posb==other.posb);
}

edge::edge(pos posa, pos posb) {
	this->posa=posa;
	this->posb=posb;
}


bool edge::is_forward() {
	return posb>posa;
}

edge_info::edge_info() {
	bp=0;
	normal_coverage=0.0;
	cancer_coverage=0.0;
	type=5;
	supporting=0;
}


//
// The state class
//

double state::score_with_flow(int flow) {
	//score regardless of edge sizes
	/*double base_unexplained=ccov-ncov;
	if (base_unexplained<0) {
		base_unexplained=-base_unexplained;
	}
	double p_unexplained=ccov-ncov*flow;
	if (p_unexplained<0) {
		p_unexplained=-p_unexplained;
	}
	double score=-(p_unexplained-base_unexplained);
	return score;*/
	//score based using edges
	/*double total=0;
	edge previous=fake_edge; edge current=fake_edge;
	for (multiset<edge>::iterator msit=gpath_set.begin(); msit!=gpath_set.end(); msit++) {
		previous=current;
		current=*msit;
		if (current==previous) {
			continue;
		} else {
			edge_info ei = edges[current];
			double ebase=abs(ei.normal_coverage-ei.cancer_coverage);
			double eflow=abs(ei.normal_coverage*flow*gpath_set.count(current)-ei.cancer_coverage);
			total+=ebase-eflow;
		}
	}*/
	//return total;
	//the funky chicken
	return scores[flow];	
}


//assumes the edge given has not been added to the state yet
void state::prepend_edge(edge e) {
	add_edge(e,false);
}

void state::append_edge(edge e) {
	add_edge(e,true);
}


void state::add_edge(edge e, bool append) {
	//fix the scoring
	add_edge_to_score(e);

	//fix up ncov,ccov,pcov
	if (gpath_set.find(e)==gpath_set.end() && gpath_set.find(e.reverse())==gpath_set.end()) {
		//add to ncov and ccov
		ncov+=edges[e].normal_coverage;
		ccov+=edges[e].cancer_coverage;	
	}	
	pcov+=edges[e].normal_coverage;

	//check if used a free edge
	if (gpath_vector.size()>0) {
		if (append) {
			edge x = gpath_vector.back();
			if (x.posb!=e.posa) {
				//used a free edge	
				fpath_set.insert(edge(x.posb,e.posa));
				if (free_edges[x.posb].find(e.posa)==free_edges[x.posb].end()) {
					cerr << "FAILED TO FIND EARLZ" << endl;
					exit(1);
				}
			}
		} else {
			edge x = gpath_vector.front();
			if (x.posa!=e.posb) {
				//used a free edge	
				fpath_set.insert(edge(e.posb,x.posa));
				if (free_edges[e.posb].find(x.posa)==free_edges[e.posb].end()) {
					cerr << "FAILED TO FIND EARLY " << e.posa.chr << ":" << e.posa.coord << " " << e.posb.chr << ":" << e.posb.coord << " " << x.posa.chr << ":" << x.posa.coord << " " << x.posb.chr << ":" << x.posb.coord << " " << endl;
					exit(1);
				}
			}
		}
	}

	//add the edge to the vector
	if (append) {
		gpath_vector.push_back(e);
	} else {
		//prepending
		gpath_vector.insert(gpath_vector.begin(),e);
	}
	//reflect the change in the set
	gpath_set.insert(e);

	//a simple cycle check DEBUG TODO
	if (gpath_vector.size()>10) {
		int sz = gpath_vector.size();
		if (gpath_vector[sz-1]==gpath_vector[sz-3] && gpath_vector[sz-3]==gpath_vector[sz-5] && gpath_vector[sz-5]==gpath_vector[sz-7]) {
			cerr << "SOMETHING WEIRD\n"; exit(1);
		}
	}


	//find the best score
	best_score();

	//make sure that this is not accidently set
	bp_check=false;
}

void state::add_edge_to_score(edge e) {
	edge_info ei = edges[e];
	double ebase=abs(ei.normal_coverage-ei.cancer_coverage);

	//find out how many times it has been used
	int count=gpath_set.count(e)+gpath_set.count(e.reverse());
	//check if the edge was already in our path
	if (count>0) {
		//remove its value from the scorings
		for (int i=MIN_FLOW; i<MAX_FLOW; i++) {
			double eflow=abs(ei.normal_coverage*i*count-ei.cancer_coverage);
			scores[i]-=ebase-eflow;
		}
	}

	//ok now lets add the thing back with a count+1
	for (int i=MIN_FLOW; i<MAX_FLOW; i++) {
		double eflow=abs(ei.normal_coverage*i*(count+1)-ei.cancer_coverage);
		scores[i]+=ebase-eflow;
	}
}

void state::best_score() {
	//find bounds
	for (set<edge>::iterator fit=fpath_set.begin(); fit!=fpath_set.end(); fit++) {
			edge fe = *fit;
			if (edges.find(fe)==edges.end()) {
				cerr << "ANOTHER BIG ERRO!" << endl;
				exit(1);
			}
			edge_info ei = edges[fe];
			/*if (free_edges[fe.posa].size()==0 || free_edges[fe.posb].size()==0) {
				cerr << "THIS IS NOT A FREE EDGE " << fe.posa.chr << ":" << fe.posa.coord << " " << fe.posb.chr << ":" << fe.posb.coord << endl;
				exit(1);
			}
			//double cancer_pairs_posa = cancer_pair_coverage[fe.posa];
			if (normal_pair_coverage.find(fe.posa)==normal_pair_coverage.end()) {
				cerr << "ANOTHER BIG ERRO! xC" << fe.posa.chr << " :" << fe.posa.coord << endl;
				exit(1);
			}
			if (normal_pair_coverage.find(fe.posb)==normal_pair_coverage.end()) {
				cerr << "ANOTHER BIG ERRO! xD" << fe.posb.chr << " :" << fe.posb.coord << endl;
				exit(1);
			}*/
			double normal_pairs_posa = normal_pair_coverage[fe.posa];
			//double cancer_pairs_posb = cancer_pair_coverage[fe.posb];
			double normal_pairs_posb = normal_pair_coverage[fe.posb];
			//int ra = 2+ceil(cancer_pairs_posa/(0.001*cancer_pairs_posa+ normal_pairs_posa));
			//int rb = 2+ceil(cancer_pairs_posb/(0.001*cancer_pairs_posb+ normal_pairs_posb));
			//max_flow=MIN(MAX(ra,rb),max_flow);
			double supporting = ((double)ei.supporting)/total_normal_pair_arcs;
			double supporting_posa = ceil(10*(supporting/(0.001*supporting + normal_pairs_posa)));
			double supporting_posb = ceil(10*(supporting/(0.001*supporting + normal_pairs_posb)));
			max_flow=MIN(MAX(supporting_posa,supporting_posb),max_flow);
	} 
	//find the best score and cp
	double best_score=0;
	int best_flow=-1;
	for (int i=MIN_FLOW; i<max_flow; i++) {
		score=score_with_flow(i);
		if (best_flow==-1 || best_score<score) {
			best_flow=i;
			best_score=score;
		}
	}
	cp=best_flow;
	score=best_score;
}


state::state() {
	//make sure some variables are set up
	scores.resize(MAX_FLOW,0.0);
	max_flow=MAX_FLOW;
	ncov=0; ccov=0; pcov=0; score=0; cp=0; 
	bp_check=false;
}

state::state(edge e) {
	scores.resize(MAX_FLOW,0.0);
	max_flow=MAX_FLOW;
	ncov=0; ccov=0; pcov=0; score=0; cp=0; 	
	
	//add in the edge
	append_edge(e);
	//find the best score
	best_score();
	//lets make sure we set this to explorable
	bp_check=false;
}


state::state(vector<edge> ve) {
	scores.resize(MAX_FLOW,0.0);
	max_flow=MAX_FLOW;
	ncov=0; ccov=0; pcov=0; score=0; cp=0; 	

	//add in all the edges
	for (unsigned int i=0; i<ve.size(); i++) {
		edge e = ve[i];
		append_edge(e);
	}

	//find the best score
	best_score();
	bp_check=false;
}

state::state(state s, edge e) {
	//lets get a copy of the parent
	*this=s; //calls copy constructor
	
	//add the edge
	append_edge(e);

	//find the best score
	best_score();

	//make sure bp_check is not set
	bp_check=false;
}


bool state::is_dup() {
	if (gpath_vector.size()==0) {
		return true;
	}


	unsigned int edges = gpath_vector.size();
	state our_state = *this; //copy our state, then we muck with it	

	unsigned int bp_so_far=0;
	edge first_e = gpath_vector[0];
	set<pos>::iterator sit;
	sit=bps.find(first_e.posa);
	while (bp_so_far<bp_range && sit->chr==first_e.posa.chr) {
		pos old_point = *sit;
		if (first_e.is_forward()) {
			sit--;
			if (sit==bps.begin()) {
				break;
			}
		} else {
			sit++;
			if (sit==bps.end()) {
				break;
			}
		}
		edge e=edge(*sit,old_point);
		bp_so_far+=e.length();
		our_state.prepend_edge(e);
	
		if (our_state.score>score) {
			//just checking
			if (gpath_vector.size()!=edges) {
				cerr << " OH BOY, error :( " << endl;
				exit(1);
			}
			return true;
		}
	}


	//just checking
	if (gpath_vector.size()!=edges) {
		cerr << " OH BOY, error :( " << endl;
		exit(1);
	}

	return false;
}


bool state::go_on() {
	return tail_check()>ZERO;
}

double state::tail_check() {
	//lets chec the last bp_range and make sure its positive
	vector<edge> ve;
	vector<edge> ve_c;
	unsigned int bp=0;

	unsigned int i=0;

	//this finds out how many edges from the back we need
	for (; i<gpath_vector.size() && bp<bp_range; i++) {
		edge e = gpath_vector[gpath_vector.size()-1-i];
		bp+=e.length();
	}
	size_t nsize = gpath_vector.size()-i;

	//make ve_c be the collection of these
	for (unsigned int n=0; n<nsize; n++) {
		ve_c.push_back(gpath_vector[n]);
	}
	//take the local path complement with respect to ve_c
	/*for (; i>0; i--) {
		edge e = gpath_vector[gpath_vector.size()-i];
		cc+=edges[e].cancer_coverage;	
		nc+=edges[e].normal_coverage;	
		ve.push_back(e);
	}*/
	
	//make the state
	//state s = state(ve);
	state s_c = state(ve_c);
	
	//cerr << "TAIL CHECK " << cc << " | " << nc << "((("  << i << " K "<< z.posa.chr << ":" << z.posa.coord << " "<< gpath_vector.size() << " vs " << nsize << " x " << this->score <<  " " << s_c.score << endl;

	//double old_score=score_with_flow(s_c.cp);	
	
	//return bp_range*(old_score-s_c.score)/length();
	return bp_range*score/length();
}

bool state::operator>(const state &other) const {
	return score>other.score;
}
bool state::operator<(const state &other) const {
	return score<other.score;
}

unsigned int state::length() {
	unsigned int l=0;
	for (unsigned int i=0; i<gpath_vector.size(); i++) {
		l+=gpath_vector[i].length();
	}
	return l;
}


string state::str() {
	ostringstream oss;
	oss << "State: " << (bp_check ? "TRUE" : "FALSE" ) << endl;
	oss << "Score: " << score << " copies: " << cp <<  " max_flow: " << max_flow << " ccov: " << ccov << " ncov: " << ncov << " length: " << length() << " fes: " << fpath_set.size() << endl;
	oss << "Tail check: " << tail_check() << " edges: " << gpath_vector.size() <<  endl;

	bool exitb=false;

	for (set<edge>::iterator fit=fpath_set.begin(); fit!=fpath_set.end(); fit++) {
		oss << "FE: " << fit->posa.chr << ":" << fit->posa.coord << "\t" << fit->posb.chr << ":" << fit->posb.coord << endl;
	}
	
	if (gpath_vector.size()>0) {
		edge start = gpath_vector[0];	
		edge last = gpath_vector[0];
		//oss << "X" << start.posa.chr << ":" << start.posa.coord << " " << start.posb.chr << ":" << start.posb.coord << endl;
		double nc = edges[start].normal_coverage;
		double cc = edges[start].cancer_coverage;
		int len = start.length();
		for (unsigned int i=1; i<gpath_vector.size(); i++) {
			edge e = gpath_vector[i];
			//if (e.posa.coord==94516290 || e.posb.coord==94516290) {
			//	exitb=true;
			//}
			//oss << "X" << e.posa.chr << ":" << e.posa.coord << " " << e.posb.chr << ":" << e.posb.coord << endl;
			if (e.posa==last.posb && (e.is_forward()==last.is_forward())) {
				//keep it going
				last=e;
				nc+=edges[e].normal_coverage;
				cc+=edges[e].cancer_coverage;
				len += e.length();
			} else {
				
				//should print and set new one
				oss << "\t" << start.posa.chr << ":" << start.posa.coord << " ~ " << last.posb.chr << ":" << last.posb.coord  << " [ " << len << " ] ";
				oss << " ncov: " << nc << " ccov: " << cc  <<  (start.posa<last.posb ? "\t+" : "\t-" ) << endl;


				edge free_edge = edge(last.posb,e.posa);
				if (edges.find(free_edge)==edges.end() || free_edges[last.posb].find(e.posa)==free_edges[last.posb].end()) {
					cerr << "Could not find the free edge!"  << endl;
					cerr << oss.str() << endl;
					exit(1);
				}
				edge_info ei = edges[edge(last.posb,e.posa)];
				if (normal_pair_coverage.find(e.posa)==normal_pair_coverage.end()) {
					cerr << "ANOTHER BIG ERRO! xCz" << e.posa.chr << " :" << e.posa.coord << endl;
					exit(1);
				}
				if (normal_pair_coverage.find(last.posb)==normal_pair_coverage.end()) {
					cerr << "ANOTHER BIG ERRO! xDz" << last.posb.chr << " :" << last.posb.coord << endl;
					exit(1);
				}
				double cancer_pairs_posa = cancer_pair_coverage[last.posb];
				double normal_pairs_posa = normal_pair_coverage[last.posb];
				double cancer_pairs_posb = cancer_pair_coverage[e.posa];
				double normal_pairs_posb = normal_pair_coverage[e.posa];
				//print some free edge info
				oss << "\ttype: " <<  ei.type << " supporting: " << ei.supporting;
				double supporting = ((double)ei.supporting)/total_normal_pair_arcs;
				double supporting_posa = ceil(10*(supporting/(0.001*supporting + normal_pairs_posa)));
				double supporting_posb = ceil(10*(supporting/(0.001*supporting + normal_pairs_posb)));
				//max_flow=MIN(3*MAX(supporting_posa,supporting_posb),max_flow);
				oss << " posa(cancer,normal) " << cancer_pairs_posa << "," << normal_pairs_posa << " RX: " << supporting_posa << " R:" << cancer_pairs_posa/(0.001*cancer_pairs_posa+ normal_pairs_posa);
				oss << " posb(cancer,normal) " << cancer_pairs_posb << "," << normal_pairs_posb << " RX: " << supporting_posb << " R:" << cancer_pairs_posb/(0.001*cancer_pairs_posb+ normal_pairs_posb) << endl;



				start=e;
				last=e;
				len=start.length();
				nc=edges[start].normal_coverage;
				cc=edges[start].cancer_coverage;
			}
		}
		edge e = gpath_vector.back();
		oss << "\t" << start.posa.chr << ":" << start.posa.coord << " ~ " << last.posb.chr << ":" << last.posb.coord  << " [ " << len << " ] ";
		oss << " ncov: " << nc << " ccov: " << cc  << (start.posa<last.posb ? "\t+" : "\t-" )  << endl;
		oss << "\t\t" << e.posa.chr << ":" << e.posa.coord << " ~ " << e.posb.chr << ":" << e.posb.coord  << " [ " << "?" << " ] ";
		oss << " ncov: " << edges[e].normal_coverage << " ccov: " << edges[e].cancer_coverage  << (start.posa<last.posb ? "\t+" : "\t-" )  << endl;
	}
	
	if (exitb) {
		cout << oss.str() << endl ;
		cout << "FOUND AND EXITED" << endl;
		exit(1);
	}

	
	return oss.str();
}


vector<state> state::children() {
	//lets move bp along the edge we are on and see where this gets us
	edge last_edge = gpath_vector.back();
	pos second_last_pos = last_edge.posa;
	pos last_pos = last_edge.posb;
	bool forward=gpath_vector.back().is_forward();
	
	
	vector<state> children;


	//get the genome path children
	if (last_edge.posa.chr!=last_edge.posb.chr) {
		cout << str() << endl;
		cerr << "MASSIVE ERROR " << endl;
		exit(1);
	}
	if (bp_check) {
		set<pos>::iterator sit = bps.find(last_pos);
		if (sit==bps.end()){
			cerr << "Something went wrong" << endl;
			exit(1);
		}	
		unsigned int bp_moved=0;

		pos current_pos=last_pos;
		if (forward) {
			sit++;
		} else {
			sit--;
		}
		while (bp_moved<bp_range && current_pos.chr==last_pos.chr && sit!=bps.end()) {
			last_pos=current_pos;
			current_pos=*sit;
			if (forward) {
				if (last_pos>current_pos) {
					cerr << "BIG ERROR X3" << endl;
					exit(1);
				}
			} else {
				if (current_pos>last_pos) {
					cerr << "BIG ERROR X4" << endl;
					exit(1);
				}

			}
			if (last_pos.chr!=current_pos.chr) {
				break;
			}
			edge e = edge(last_pos,current_pos);
			
			state child;
			//want to add child with this edge on it
			if (children.size()==0) {
				child=state(*this,e);
			} else {
				child=state(children.back(),e);
			}
			children.push_back(child);
		
			if (forward) {
				sit++;
			} else {
				sit--;
			}
			
			bp_moved+=e.length();
		}
		if (children.size()>0) {
			children.back().bp_check=true;
		}
	}
	
	//in any case lets check for donor edges we can add from the current state!
	last_edge = gpath_vector.back();
	second_last_pos = last_edge.posa;
	last_pos = last_edge.posb;
	for (set<pos>::iterator sit = free_edges[last_pos].begin(); sit!=free_edges[last_pos].end(); sit++) {
		edge fe = edge(last_pos,*sit);
		if (fpath_set.find(fe)!=fpath_set.end() || fpath_set.find(fe.reverse())!=fpath_set.end() || fpath_set.size()==MAX_FREE) {
			//already used this free edge
			continue;
		}
		edge_info fei = edges[fe];
		//can enter the edge?
		//
		int type = fei.type;
		int supporting = fei.supporting;
		if (supporting==0) {
			cerr << "HUGE ERROR"  << fe.posa.chr << ":" << fe.posa.coord << " " << fe.posb.chr << ":" << fe.posb.coord << endl;	
			exit(1);
		}
		set<pos>::iterator ssit;

		edge to_add = fake_edge;

		//check to see how we can use the edge
		if (forward && type%2==0) {
			if (type==0) {
				//type 0 - leave on positive
				ssit=bps.find(*sit);
				ssit++;
				if (ssit!=bps.end() && sit->chr==ssit->chr) {
					//use this edge
					to_add = edge(*sit,*ssit);
				}
			} else {
				//type 2 - leave on negative
				ssit=bps.find(*sit);
				if (ssit!=bps.begin()) {
					ssit--;
					//use this edge
					if (sit->chr==ssit->chr) {
						to_add = edge(*sit,*ssit);
					}
				}
			}
		} else if (!forward && type%2==1) {
			if (type==1) {
				//type 1 - leave on negative
				ssit=bps.find(*sit);
				if (ssit!=bps.begin()) {
					ssit--;
					//use this edge
					if (sit->chr==ssit->chr) {
						to_add = edge(*sit,*ssit);
					}
				}
			} else {
				//type 3 - leave on positive
				ssit=bps.find(*sit);
				ssit++;
				if (ssit!=bps.end() && sit->chr==ssit->chr) {
					//use this edge
					to_add = edge(*sit,*ssit);
				}
			} 

		}

		//check if we found a edge to use, if so add it!	
		if (!(to_add==fake_edge)) {
			state child = state(*this,to_add);
			child.bp_check=true;
			children.push_back(child);
		}
	}

	return children;	
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

//read in the bp arcs
void read_arcs(char * filename, bool normal) {
	ifstream f (filename);
	string chr_s;
	unsigned int coord;
	int arcs;
	unsigned int total=0;
	while (f) {
		string line;
		getline(f,line);
		istringstream is(line);
		is >> chr_s >> coord >> arcs;
		int chr = to_chr(chr_s.c_str());
		
		pos p = pos(chr,coord);
		if (normal) {		 
			normal_pair_coverage[p]=arcs;
		} else {
			cancer_pair_coverage[p]=arcs;
		}
		total+=arcs;
	}

	set<pos>::iterator sit;
	if (normal) {
		for (map<pos,double>::iterator sit = normal_pair_coverage.begin(); sit!=normal_pair_coverage.end(); sit++) {
			normal_pair_coverage[sit->first]=sit->second/total;
		}
		total_normal_pair_arcs=total;
	} else {
		for (map<pos,double>::iterator sit = cancer_pair_coverage.begin(); sit!=cancer_pair_coverage.end(); sit++) {
			cancer_pair_coverage[sit->first]=sit->second/total;
		}
		total_cancer_pair_arcs=total;
	}
	

	cerr << "Read " << total << " arcs from " << filename << endl;

}

//read in the edges from clustering
void read_links(char * filename) {
	ifstream f (filename);
	//char chr[10]="";
	unsigned int bpa,bpb,l_from,l_to,cluster_idx;
	double avg_md;
	string chra,chrb;
	int type,total;
	int total_links=0;
	while (f) {
		string line;
		getline(f,line);
		istringstream is(line);
		//chr2    144365053       185125297       3       2       2       0       2.8396e+08      1012    chr3    EDGE
		//chr3    84440059        190840923       3       9       67      70      3.04423e+08     1013    chr4    EDGE
		is >> chra >> bpa >> bpb >> type >> total >> l_from >> l_to >> avg_md >> cluster_idx >> chrb; 
		//cout << chra << bpa << bpb << type << total << l_from << l_to << avg_md << chrb << endl;
		//cout << line << endl << line << endl;
		int nchra=to_chr(chra.c_str());
		int nchrb=to_chr(chrb.c_str());
		if (bpa<=0 || bpb<=0) {
			continue;
		}
		total_links++;
		pos posa=pos(nchra,bpa);
		pos posb=pos(nchrb,bpb);

		bps.insert(posa);
		bps.insert(posb);

		//add the one direction
		edge ea = edge(posa,posb);

		free_edges[posa].insert(posb);
		edges[ea].type=type;
		edges[ea].supporting=total;
		edges[ea].bp=ea.length();

		//add the other direction
		edge eb = ea.reverse();

		if (type<=1) {
			type=1-type;
		}
		free_edges[posb].insert(posa);
		edges[eb].type=type;
		edges[eb].supporting=total;
		edges[eb].bp=eb.length();
	}
	cerr << "Read " << total_links << " links from " << filename << endl;	
	return ;
}




void read_cov(char * filename, bool normal) {
	cerr << "Reading coverage from file " << filename << endl;
	gzFile gzf = gzopen(filename,"r");
	if (gzf==NULL) {
		fprintf(stderr, "Failed to open file %s\n",filename);
		exit(1);
	}
	
	//find the size of one entry
	size_t soe = sizeof(unsigned short)+sizeof(unsigned int)+sizeof(unsigned short);
	size_t chunk = 1024*1024*1024;
	size_t size_so_far = 0;
	char * buffer = (char*) malloc(chunk);
	if (buffer==NULL) {
		cerr << " FALLED TO MALLOC " << endl;
		exit(1);
	}

	//the read shunt
	/*for (int i=0; i<6; i++) {
		int read = gzread(gzf,buffer+size_so_far,chunk);
		size_so_far+=read;
		cerr << "Warning!!!!" << endl;
		buffer=(char*)realloc(buffer,size_so_far+chunk);
	}
	size_so_far=(size_so_far/soe)*soe;*/

	//the real read loop
	while (!gzeof(gzf)) {
		int read = gzread(gzf,buffer+size_so_far,chunk);
		cerr << "\rRead so far " << size_so_far;
		if (read<0) {
			cerr << "Error reading gzipped file!\n";
			exit(1);
		}
		size_so_far+=read;
		buffer=(char*)realloc(buffer,size_so_far+chunk);
		if (buffer==NULL) {
			cerr << " FALLED TO REALLOC " << endl;
			exit(1);
		}
	}
	
	//lets get a buffer to fit the file	
	cerr << "Done reading file  " << size_so_far <<  endl;

	size_t sz = size_so_far;	

	unsigned int entries = sz/soe;
	map<pair<int,unsigned int> , double> m;
	unsigned long total_coverage=0;


	map<edge, edge_info > edges_ret;
	
	#pragma omp parallel 
	{
	unsigned int threads = omp_get_num_threads();
	unsigned int thread_id = omp_get_thread_num();
	//cerr << "thread " <<  thread_id << endl;
	unsigned long total_coverage_t=0;

	map<edge, edge_info > edges_t;

	set<pos>::iterator it = bps.begin();
	pos prev = *it; it++;

	unsigned short chr, cov;
	unsigned int coord;
	for (unsigned int i=0; i<entries; i++) {
		if (i%threads!=thread_id) {
			continue; // not our job!
		}
		char* base = buffer+i*soe;
		chr=*((unsigned short *)base);
		base+=sizeof(unsigned short);
		coord=*((unsigned int *)base);
		base+=sizeof(unsigned int);
		cov=*((unsigned short *)base);

		total_coverage_t+=cov;

		
		pos p = pos(chr,coord);
		while (p>*it && it!=bps.end()) {
			prev=*it;
			it++;
		}
		if (it==bps.end()) {
			//cout << "broke at " << p.first << " " << p.second << "    " << (*it).first << " " << (*it).second << endl;
			break;	
		}
		if (prev.chr==chr) {
			edge ea = edge(prev,*it);
			edge eb = ea.reverse();
			if (normal) {
				edges_t[ea].normal_coverage+=cov;
				edges_t[eb].normal_coverage+=cov;
			} else {
				edges_t[ea].cancer_coverage+=cov;
				edges_t[eb].cancer_coverage+=cov;
			}
		}

	}

	#pragma omp critical 
	{
		//iterator over thread specific and add to main
		total_coverage+=total_coverage_t;
		for (map<edge,edge_info>::iterator mit=edges_t.begin(); mit!=edges_t.end(); mit++) {
			edges[mit->first].normal_coverage+=mit->second.normal_coverage;
			edges[mit->first].cancer_coverage+=mit->second.cancer_coverage;
		}
	}
	} //end openmp section

	for (map<edge,edge_info>::iterator mit=edges.begin(); mit!=edges.end(); mit++) {
		if (normal) {
			edges[mit->first].normal_coverage/=total_coverage;
		} else {
			edges[mit->first].cancer_coverage/=total_coverage;
		}
	}

	cerr << "total: " << total_coverage << endl;
	free(buffer);
}

int main(int argc, char ** argv) {
	//need to load in files
	if (argc!=6) {
		printf("%s links cov_cancer cov_normal pairs_cancer pairs_normal\n", argv[0]);
		exit(1);
	}

	char * links_filename=argv[1];
	char * cov_cancer_filename=argv[2];
	char * cov_normal_filename=argv[3];
	char * pairs_cancer_filename=argv[4];
	char * pairs_normal_filename=argv[5];
	omp_set_num_threads(THREADS);
	//read in the free edges
	read_links(links_filename);

	fake_edge=edge(pos(0,0),pos(0,0));

	//read in the arc weights
	read_arcs(pairs_normal_filename,true);	
	read_arcs(pairs_cancer_filename,false);	

	//lets slice up the rest
	set<pos> to_add;
	for (set<pos>::iterator it=bps.begin(); it!=bps.end(); ) {
		pos current=*it;
		//cout << current.chr << ":" << current.coord << endl;
		it++;
		pos next=*it;
		if (it==bps.end()) {
			break;
		}
		if (current.chr==next.chr && next.coord-current.coord>MAX_EDGE_SIZE) {
			pos i = pos(current.chr,current.coord+MAX_EDGE_SIZE);
			while (next.coord>MAX_EDGE_SIZE+i.coord) {
				to_add.insert(i);
				i = pos(current.chr,i.coord+MAX_EDGE_SIZE);
			}
		}
	}
	for (set<pos>::iterator it=to_add.begin(); it!=to_add.end(); it++) {
		bps.insert(*it);
	}

	
	//initialize the rest of the edges, so can multithread
	for (set<pos>::iterator it=bps.begin(); it!=bps.end(); ) {
		free_edges[*it].size();		

		pos current=*it;
		//cout << current.chr << ":" << current.coord << endl;
		it++;
		pos next=*it;
		if (it==bps.end()) {
			break;
		}
	
		//else have current and next , lets figure it out
		if (current.chr==next.chr) {
			edge ea = edge(current,next);
			edges[ea].bp=ea.length();
			edge eb = ea.reverse();
			edges[eb].bp=eb.length();
		}
	}

	//read in the normal
	read_cov(cov_normal_filename,true);
	read_cov(cov_cancer_filename,false);

	//read in the pairs coverage
	cerr << "Done reading in data..." << endl;

	set<pos>::iterator sit = bps.begin(); 
	pos previous;
	pos current=*sit;
	sit++;
	vector<edge> start_edges;
	while (sit!=bps.end()) {
		previous=current;
		current=*sit;
		if (current.chr==previous.chr) {
			edge ea = edge(previous,current);
			edge eb = ea.reverse();
			if (free_edges[ea.posa].size()>0 || free_edges[ea.posb].size()>0) {
				start_edges.push_back(ea);
				start_edges.push_back(eb);
			}
			//cerr << "Warning only forward edges " << endl;
		}
		sit++;
	}

	cerr << "Starting enumeration..." << endl;

	//for debuggin
	/*for (int i=0; i<start_edges.size(); i++) {
		edge e = start_edges[i];
		if (e.posa.coord!=201862800 || !e.is_forward()) {
			continue;
		}
		state s = state(e);
		set<pos>::iterator bps_it = bps.find(e.posb);
		bps_it++; bps_it++;
		while (s.length()<9000000 && bps_it!=bps.end() && (*bps_it).chr==1) {
			cout << s.str();
			s = state(s,edge(s.gpath_vector.back().posb,*bps_it));
			bps_it++;
		}
	}
	exit(1);*/	

	//The main loop
	state best_state;
	set<edge> skip_edges;	


	int iteration = 0;
	do {
		iteration++;

		cerr << endl << iteration << "-th iteration , skipping : " << skip_edges.size()  << endl;
		best_state=state();
		#pragma omp parallel for schedule(dynamic,1)
		for (unsigned int i=0; i<start_edges.size(); i++) {
			//if (e.posa.coord!=201862800 || !e.is_forward()) {
			//	continue;
			//}
			edge e = start_edges[i];
			bool skip=false;
			#pragma omp critical
			{
				if (skip_edges.find(e)!=skip_edges.end()) {
					skip=true;
				}		
				cerr << "\r" << i << " thread " << omp_get_thread_num() << "   ";
			}
			if (skip) {
				continue;
			}

			priority_queue<state> pq;
			state s = state(e);
			s.bp_check=true;
			map<int, state> best_states;
			if (!s.is_dup()) { 
				best_states[0]=s;
			}
			//cout << s.str() <<  s.score << endl;
			if (s.score>ZERO && s.go_on()) {
				pq.push(s);
			}


			set<state_hash> visited;

			while (!pq.empty()) {
				state s = pq.top();
				pq.pop();

				state_hash sh = state_hash(s);

				if (visited.count(sh)>0) {
					continue;
				} else {
					//mark this as visited : TODO right way to handle?
					visited.insert(sh);
				}
				
				//cout << "CHILDREN OF " << endl << s.str() << endl;

				if (!s.is_dup() && s.score>best_states[s.fpath_set.size()].score) {
					best_states[s.fpath_set.size()]=s;
				}

				vector<state> children = s.children();

				for (unsigned int i=0; i<children.size(); i++) {
					state c = children[i];
					if (c.score>ZERO && c.go_on()) {
						//cout << "CHILD" << endl << c.str() << endl;
						pq.push(c);
					} else {
						//cout << "CHILD XXX" << endl << c.str() << endl;
					}		
				}	
			}
		

			#pragma omp critical
			{
				for (map<int, state>::iterator mit=best_states.begin(); mit!=best_states.end(); mit++) {
					if ((mit->second).score>ZERO) {
						//cout << (mit->second).str();
						if ((mit->second).score>best_state.score) {
							best_state=mit->second;			
						}
						if ((mit->second).score>edge_bests[e].score) {
							edge_bests[e]=(mit->second);
						}
					}
				}
			}
		}

		//output the best path
		cout << best_state.str() << endl;
	

		//lets change out the path
		for (unsigned int i=0; i<best_state.gpath_vector.size(); i++) {
			edge e = best_state.gpath_vector[i];
			edge_info ei = edges[e];
			ei.cancer_coverage -= ei.normal_coverage*best_state.cp;	
			edges[e]=ei;
		
			edge er = e.reverse();
			edge_info eri = edges[e];
			eri.cancer_coverage -= eri.normal_coverage*best_state.cp;			
			edges[er]=eri;
		}

		//and remove the free edges
		for (multiset<edge>::iterator msi=best_state.fpath_set.begin(); msi!=best_state.fpath_set.end(); msi++) {
			edge e = *msi;
			free_edges[e.posa].erase(e.posb);
			free_edges[e.posb].erase(e.posa);
		}

			
		//lest find out which ones we need to rerun	
		skip_edges.empty();
		#pragma omp parallel for schedule(dynamic,1)
		for (unsigned int i=0; i<start_edges.size(); i++) {
			edge e = start_edges[i];
			state edge_best;
			#pragma omp critical
			{
				edge_best = edge_bests[e];
			}
			//check if any edges overlap
			bool overlap=false;
			for (multiset<edge>::iterator sit=best_state.gpath_set.begin(); sit!=best_state.gpath_set.end(); sit++) {
				edge ea = *sit;
				edge eb = e.reverse();
				for (multiset<edge>::iterator ssit=edge_best.gpath_set.begin(); ssit!=edge_best.gpath_set.end(); ssit++) {
					if (*ssit==ea || *ssit==eb) {
						overlap=true;
					}
				}
			}
			if (overlap) {

			} else {
				#pragma omp critical
				{
					skip_edges.insert(e);	
				}
			}	
		}
	} while (best_state.score>ZERO);
	
	return 0;
}
