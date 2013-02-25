#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<set>
#include <fstream>
#include <sstream>
#include <map>
#include <string.h>
#include <vector>
#include <queue>
#include <omp.h>


#define MAX_ITERATIONS	1000000
#define MIN_FLOW 4
#define MAX_FREE 15
#define THREADS 20
#define MIN_MPAIR	4

using namespace std;

map<pair<int, unsigned int> , double> n_cov;
map<pair<int, unsigned int> , double> c_cov;
set< pair<int, unsigned int> > bps;
map< pair<int, unsigned int> , set< pair<int, unsigned int> > > edges;
map< pair<int, unsigned int> , set< pair<int, unsigned int> > > sv_edges;
map< pair< pair<int, unsigned int> , pair<int, unsigned int > > , pair<int,int> > sv_edge_types;

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


void read_links(char * filename) {
	ifstream f (filename);
	//char chr[10]="";
	unsigned int bpa,bpb,l_from,l_to,cluster_idx;
	double avg_md;
	string chra,chrb;
	int type,total;
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
		pair<int,unsigned int> pair_bpa = pair<int,unsigned int>(nchra,bpa);
		pair<int,unsigned int> pair_bpb = pair<int,unsigned int>(nchrb,bpb);
		bps.insert(pair_bpa);
		bps.insert(pair_bpb);
		sv_edges[pair_bpa].insert(pair_bpb);
		sv_edges[pair_bpb].insert(pair_bpa);
		sv_edge_types[pair< pair<int, unsigned int> , pair<int, unsigned int > >(pair_bpa,pair_bpb)]=pair<int,int>(type,total);
		if (type<=1) {
			sv_edge_types[pair< pair<int, unsigned int> , pair<int, unsigned int > >(pair_bpb,pair_bpa)]=pair<int,int>(1-type,total);
		} else {
			sv_edge_types[pair< pair<int, unsigned int> , pair<int, unsigned int > >(pair_bpb,pair_bpa)]=pair<int,int>(type,total);
		}
	}
	return ;
}

void usage(char * p) {
	fprintf(stderr,"%s links cov_normal cov_control\n",p);
}

map<pair<int, unsigned int> , double> read_cov(char * filename) {
	FILE * fptr = fopen(filename,"rb");
	if (fptr==NULL) {
		fprintf(stderr, "An error opening %s has occured\n",filename);
	}


	//get the file size
	size_t soe = sizeof(unsigned short)+sizeof(unsigned int)+sizeof(unsigned short);
	cerr << " Start reading file  " << endl;
	fseek(fptr, 0L, SEEK_END);
	size_t sz = ftell(fptr);
	if (sz%soe!=0) {
		cerr << "FILE CORRUPT!\n";
		exit(1);
	}
	fseek(fptr, 0L, SEEK_SET);
	char * buffer = (char*) malloc(sz);
	if (buffer==NULL) {
		cerr << " FALLED TO MALLOC " << endl;
		exit(1);
	}
	fread(buffer,sz,1,fptr);
	cerr << " Done reading file  " << endl;

	unsigned int entries = sz/soe;
	map<pair<int,unsigned int> , double> m;
	unsigned long total_coverage=0;
	omp_set_num_threads(THREADS);
	#pragma omp parallel 
	{
	unsigned short chr, cov;
	unsigned int pos;
	pair<int, unsigned int> p;
	set<pair<int, unsigned int> >::iterator it = bps.begin();
	pair<int, unsigned int> prev=*it; it++;
	int threads = omp_get_num_threads();
	int thread_id = omp_get_thread_num();
	cerr << "thread " <<  thread_id << endl;
	//need to merge these after
	map<pair<int,unsigned int> , double> m_t;
	unsigned long total_coverage_t=0;
	for (unsigned int i=0; i<entries; i++) {
		if (i%threads!=thread_id) {
			continue; // not our job!
		}
		char* base = buffer+i*soe;
		chr=*((unsigned short *)base);
		base+=sizeof(unsigned short);
		pos=*((unsigned int *)base);
		base+=sizeof(unsigned int);
		cov=*((unsigned short *)base);

		total_coverage_t+=cov;
		p.first=chr; p.second=pos;
		while (p>*it && it!=bps.end()) {
			prev=*it;
			it++;
		}
		if (it==bps.end()) {
			//cout << "broke at " << p.first << " " << p.second << "    " << (*it).first << " " << (*it).second << endl;
			break;	
		}
		if (prev.first==chr) {
			m_t[prev]+=cov;
		}

	}

	#pragma omp critical
	{
		for (map<pair<int,unsigned int> , double>::iterator it = m_t.begin(); it!=m_t.end(); it++) {
			m[it->first]+=it->second;
		}
		total_coverage+=total_coverage_t;
	}
	} //end openmp section

	free(buffer);
	/*while(!feof(fptr)) {
		fread(&chr,sizeof(unsigned short),1,fptr);
		fread(&pos,sizeof(unsigned int),1,fptr);
		fread(&cov,sizeof(unsigned short),1,fptr);
		//cout << chr << " " << pos << " C:" << cov << endl;
	}*/

	//lets normalize
	map< pair<int, unsigned int> , double>::iterator itm = m.begin();
	for (; itm!=m.end(); itm++) {
		m[itm->first]=m[itm->first]/total_coverage;
		//cout << m[itm->first] << endl;
	}


	return m;
}


class State {
	public:
		vector< pair< pair<int, unsigned int>, pair<int, unsigned int> > > edges_used;
		double score;
		double copy_count;
		int free_edges;
		State();
		State(State, pair< pair<int, unsigned int>, pair<int, unsigned int> >); //add this edge
		State extend_normal();
		void print();
		void best_score();
		bool is_cycle(bool p);
		bool operator<(const State &other) const;
		bool operator>(const State &other) const;
};

State::State() {
	score=-10000;
	copy_count=0;
	free_edges=0;
	return;
}

State::State(State s, pair<pair<int, unsigned int> , pair<int,unsigned int> > e) {
	edges_used=vector< pair< pair<int, unsigned int>,pair<int, unsigned int> > >(s.edges_used);
	free_edges=s.free_edges;
	//find out if this is a free edge
	if (edges_used.size()>0) {
		pair<int, unsigned int> last_edge = edges_used.back().second;
		if (e.first!=last_edge) {
			//its a free edge
			//find out if it has been used before
			pair<pair<int, unsigned int> , pair<int, unsigned int> > free_e;
			free_e.first=last_edge;
			free_e.second=e.first;
			pair<pair<int, unsigned int> , pair<int, unsigned int> > free_e_inv;
			free_e_inv.first=free_e.second;
			free_e_inv.second=free_e.first;
			for (int i=1; i<edges_used.size(); i++){ 
				pair<pair<int, unsigned int> , pair<int, unsigned int> > free_e_used;
				free_e_used.first=edges_used[i-1].second;
				free_e_used.second=edges_used[i].first;
				if (free_e_used==free_e || free_e_used==free_e_inv) {
					//edge has been already used, lets drop this
					return;
				}		
			}
			
		}
	}
	edges_used.push_back(e);
	best_score();//set the best score
	return;
}

bool State::operator<(const State &other) const {
	return score<other.score;
}

bool State::operator>(const State &other) const {
	return score>other.score;
}



void State::print() {
	cout << "score: " << score << "," << copy_count << " ";
	cout << "Edges: " ;
	vector <pair<int, unsigned int> > v;
	for (int i=0; i<edges_used.size(); i++) {
		pair<int, unsigned int> ea=edges_used[i].first;
		pair<int, unsigned int> eb=edges_used[i].second;
		unsigned int s=v.size();
		if (s==0) {
			v.push_back(ea);
			v.push_back(eb);
		} else {
			if (v.back()==ea) {
				if (ea==v[s-1] && eb==v[s-2]) {
					cout << "ERROR2RT$"  << endl;
					exit(1);
				}
				v[s-1]=eb;
			} else {
				v.push_back(ea);
				v.push_back(eb);
			}
		}
	}
	bool range=true;
	for (int i=0; i<v.size(); i++) {
		pair<int, unsigned int> p = v[i];
		int l=0;
		if (i%2==1) {
			pair<int , unsigned int> n = v[i-1];
			l=abs(n.second-p.second);
			cout << p.first << ":" << p.second << "[" << l << "]" << (range ? "~" : "<->");
		} else {
			cout << p.first << ":" << p.second << (range ? "~" : "<->");
		}
		range=!range;
	}
	/*cout << " |||  "; 
	for (int i=0; i<edges_used.size(); i++) {
		pair<int, unsigned int> bp1=edges_used[i].first;
		pair<int, unsigned int> bp2=edges_used[i].second;
		pair<int ,unsigned int> first=bp1;
		if (bp2<bp1) {
			first=bp2;
		}
		cout << bp1.first << ":" << bp1.second << "-" << bp2.first << ":" << bp2.second << " c: " << c_cov[first] << " n: " << n_cov[first] << ",";
	}*/
	cout << endl;
}

State State::extend_normal() {
		pair<int, unsigned int> last_bp = edges_used.back().second;
		pair<int, unsigned int> second_last_bp = edges_used.back().first;
		pair< pair<int, unsigned int> , pair<int, unsigned int> > e;
		set<pair<int, unsigned int> >::iterator fit;
		//add in the normal edge
		if (last_bp>second_last_bp) {
			//add forward edge
			fit=bps.find(last_bp);
			fit++;
			if (fit!=bps.end()) {
				pair<int, unsigned int> next_bp=*fit;
				e.first=last_bp;
				e.second=next_bp;
				State ns=State(*this,e);
				return ns;
			}
		} else {
			//add backward edge
			fit=bps.find(last_bp);
			if (fit!=bps.begin()) {
				fit--;
				pair<int, unsigned int> prev_bp=*fit;
				e.first=last_bp;
				e.second=prev_bp;
				State ps=State(*this,e);
				return ps;
			}
		}
		return State();
}


bool State::is_cycle(bool p) {
	pair<int, unsigned int> first_bp = edges_used.front().first;
	pair<int, unsigned int> last_bp = edges_used.back().second;
	pair<int, unsigned int> second_last_bp = edges_used.back().first;
	//iterate over all sv_edges
	if (sv_edges[last_bp].find(first_bp)!=sv_edges[last_bp].end()) {
		//have found it!
		return true;
	}
	//check for existence 
	if (first_bp==last_bp) {
		return true;
	}	
	/*set<pair<int, unsigned int> >::iterator fit;
	if (last_bp>second_last_bp) {
		//add forward edge
		fit=bps.find(last_bp);
		fit++;
		if (*fit==first_bp) {
			if (p) {
				cout << "X2";
			}
			return true;
		}
	} else {
		//add backward edge
		fit=bps.find(last_bp);
		if (fit!=bps.begin()) {
			fit--;
			if (*fit==first_bp) {
				if (p) {
					cout << "X4";
				}
				return true;
			}
		}
	}
	if (p) {
		cout << "X3";
	}*/
	return false;
		
}

void State::best_score() {
	//try between 1 and 50
	//find out how many times each free edge was used and if we have enough to use them
	int local_max_flow=100;
	map<pair<pair<int, unsigned int>,pair<int, unsigned int> > , int> z;
	for (int i=1; i<edges_used.size(); i++) {
		pair<pair<int, unsigned int>,pair<int, unsigned int> > current = edges_used[i];
		pair<pair<int, unsigned int>,pair<int, unsigned int> > prev = edges_used[i-1];
		if (prev.second.first!=current.first.first || prev.second.second!=current.first.second) {
			//its a free edge
			pair<pair<int, unsigned int>,pair<int, unsigned int> > k;
			if (prev.second<current.first) {
				k=pair<pair<int, unsigned int>,pair<int, unsigned int> >(prev.second,current.first);
			} else {
				k=pair<pair<int, unsigned int>,pair<int, unsigned int> >(current.first,prev.second);
			}
			int total=sv_edge_types[k].second;
			z[k]+=1;
			int mx = total/(z[k]*MIN_MPAIR);
			if (mx<local_max_flow) {
				local_max_flow=mx;
			}
			
		}	
	}
	//find out how many times each edge is used
	map<pair<int, unsigned int>, int > ms;
	int max_edge_count=0;
	double n_total=0.0;
	double c_total=0.0;
	double npath_total=0.0;
	int ctx=0;
	int last_chr=0;
	for (vector<pair<pair<int, unsigned int>,pair<int, unsigned int> > >::iterator it=edges_used.begin(); it!=edges_used.end(); it++) {
		pair<int, unsigned int> pa=(*it).first;
		pair<int, unsigned int> pb=(*it).second;
		if (last_chr==0) {
			last_chr=pb.first;
		} else {
			if (last_chr!=pa.first) {
				ctx++;
				last_chr=pb.first;
			}
			last_chr=pb.first; 
		}
		if (pa.first==pb.first) {
			if (pa>pb) {
				if (ms[pb]==0) {
					c_total+=c_cov[pb];
					n_total+=n_cov[pb];
				}
				npath_total+=n_cov[pb];
				ms[pb]+=1;
				if (ms[pb]>max_edge_count) {
					max_edge_count=ms[pb];
				}
			} else {
				if (ms[pa]==0) {
					c_total+=c_cov[pa];
					n_total+=n_cov[pa];
				}
				npath_total+=n_cov[pa];
				ms[pa]+=1;
				if (ms[pa]>max_edge_count) {
					max_edge_count=ms[pa];
				}
			}
		} else {
			//sv edge
		}
	}
	//compute the function
	vector<double> v (100,0.0);
	for (int i=0; i<100; i++) {
		if (i>local_max_flow || (i!=0 && i<MIN_FLOW)) {
			v[i]=0;
		} else {
			double i_unexplained=c_total-npath_total*i;
			if (i_unexplained<0) {
				i_unexplained=-i_unexplained;
			}
			double base_unexplained=c_total-n_total;
			if (base_unexplained<0) {
				base_unexplained=-base_unexplained;
			}
			v[i]=v[i]-(i_unexplained-base_unexplained);
		}
	}
	/*for (map<pair<int,unsigned int> , int >::iterator it=ms.begin(); it!=ms.end(); it++) {
		//cout << it->first.first << " " << it->first.second << " " << it->second << endl;
		if (c_cov.find(it->first)==c_cov.end()) {
			cout << "BIG ERROR!" << it->first.first << ":" << it->first.second << endl;
		}
		if (n_cov.find(it->first)==n_cov.end()) {
			cout << "BIG ERROR!" << it->first.first << ":" << it->first.second << endl;
		}
		double c = c_cov[it->first];
	        double n = n_cov[it->first];	
		for (int i=0; i<100; i++) {
			if (i!=0 && i<MIN_FLOW) {
				v[i]=0;
				continue;
			}
			double i_unexplained=c-n*i*it->second;
			if (i_unexplained<0) {
				i_unexplained=-i_unexplained;
			}
			double base_unexplained=c-n;
			if (base_unexplained<0) {
				base_unexplained=-base_unexplained;
			}
			v[i]=v[i]-(i_unexplained-base_unexplained);
			//cout << i << " " << v[i] << endl;
		}
	}*/
	int best_i=0;
	double best_v=v[0];
	for (int i=0; i<100; i++) {
		if (v[i]>best_v) {
			best_v=v[i];
			best_i=i;
		}
	}

	//pretty much forced a min flow 
	if (best_i==MIN_FLOW) {
		score=0;
		copy_count=0;
		return;
	}

	for (int i=0; i<free_edges; i++) {
		score=score*0.95;
	}
	for (int i=0; i<ctx; i++) {
		score=score*0.9;
	}
	
	//set the score
	score=best_v;
	copy_count=best_i;
}

int main(int argc, char ** argv) {
	if (argc!=4) {
		usage(argv[0]);
		exit(1);
	}


	omp_set_num_threads(THREADS);

	char * links_filename=argv[1];
	char * normal_cov_filename=argv[2];
	char * cancer_cov_filename=argv[3];

	//have the sv edges
	read_links(links_filename);


	//have regular edges
	set< pair<int, unsigned int> >::iterator it=bps.begin();
	pair<int, unsigned int> prev=*it; it++;
	for (; it!=bps.end(); it++) {
		pair<int, unsigned int> p = *it;
		int chr=p.first;
		unsigned int pos=p.second;
		if (prev.first==chr) {
			//add an edge from previous to here
			edges[prev].insert(p);
			edges[p].insert(prev);
		}	
		prev=*it;
		//cout << chr << " " << pos << endl;
	}

	//read in normal coverage
	cout << "loading coverage..." << endl;
			n_cov = read_cov(normal_cov_filename);
			c_cov = read_cov(cancer_cov_filename);
	/*
	#pragma omp parallel for
	for (int i=0; i<2; i++) {
		if (i==0) {
			n_cov = read_cov(normal_cov_filename);
		} else if (i==1) {
			c_cov = read_cov(cancer_cov_filename);
		}
	}*/


	double sum=0;
	for (set<pair<int, unsigned int> >::iterator sit = bps.begin(); sit!=bps.end(); sit++) {
		sum+=c_cov[*sit];
		sum+=n_cov[*sit];
	}

	cout << "Begun search..." << endl;
	for (set<pos>::iterator zi=free_edges[pos(12,41183490)].begin(); zi!=free_edges[pos(12,41183490)].end(); zi++) {
		cout << zi->chr << " : " << zi->coord << endl;
	}

	exit(1);
	vector<pair<int,unsigned int> > v (bps.begin(),bps.end());
	#pragma omp parallel for schedule(dynamic, 1)
	for (int i=0; i<v.size(); i++) {
		int thread_id = omp_get_thread_num();
		//cerr << "S:" << thread_id << ":" << i << endl;
		//itialize the state
		State o = State();
		State best_state = o;
		State best_cycle = o;
		//initialize the q
		priority_queue<State> pq;
		pair< pair<int, unsigned int> , pair<int, unsigned int> > e;
		for (set<pair<int,unsigned int> >::iterator sit=edges[v[i]].begin(); sit!=edges[v[i]].end(); sit++) {
			pair<int, unsigned int> other = *sit;
			if (other.first==v[i].first) {
				//same chr good news its an edge!
				e.first=v[i];
				e.second=other;
				State x = State(o,e);		
				pq.push(x);
				if (x.score>best_state.score) {
					best_state=x;
				}
				if (x.is_cycle(false) && x.score>best_cycle.score) {
					best_cycle=x;
				}
			}
		}


		int iterations=0;
		//lets pop the queue
		while (!pq.empty()) {
			if (iterations>MAX_ITERATIONS) {
				cout << "REACHED MAX IT" << endl;
				break;
			}
			//cerr << thread_id << ":looping:" << i << endl;
			iterations++;
			State x=pq.top();
			pq.pop();
			if (x.score<=10e-14 || x.free_edges>MAX_FREE) {
				continue;
			}
			if (x.score>best_state.score) {
				best_state=x;
			}
			if (x.is_cycle(false) && x.score>best_cycle.score) {
				best_cycle=x;
			}
			//cout << " score: " << x.score << " cc: " << x.copy_count << endl;
			//x.print();
			//lest seed the next generation
			pair<int, unsigned int> last_bp = x.edges_used.back().second;
			pair<int, unsigned int> second_last_bp = x.edges_used.back().first;
			bool forward=true;
			if (last_bp.second<second_last_bp.second) {
				forward=false;
			}
			pair< pair<int, unsigned int> , pair<int, unsigned int> > e;
			set<pair<int, unsigned int> >::iterator fit;

			pair<pair<int, unsigned int>, pair<int, unsigned int> > first_edges=x.edges_used.front();
			//first lets use free edges
			for (set< pair<int, unsigned int> >::iterator it=sv_edges[last_bp].begin(); it!=sv_edges[last_bp].end(); it++) {
				//lets add forward and backward edges
				pair<int, unsigned int> other_bp = *it;
				if (sv_edge_types.find(pair< pair<int, unsigned int> , pair<int, unsigned int > >(last_bp,other_bp))==sv_edge_types.end()) {
					cerr << "MISSING FREE TYPE!\n";
					exit(1);
				}
				pair<int,int> pi=sv_edge_types[pair< pair<int, unsigned int> , pair<int, unsigned int > >(last_bp,other_bp)];
				int type=pi.first;
				int support=pi.second;
				if (forward && type%2==1) {
					continue;
				} else if (!forward && type%2==0) {
					continue;
				}
				if (type>=2) {
					forward=!forward;
				}
				//lets get the next and previous
				fit = bps.find(other_bp);
				fit++;
				//TODO CHECK FOR SAME CHR!!
				if (fit!=bps.end() && forward) {
					pair<int, unsigned int> next_bp=*fit;
					e.first=other_bp;
					e.second=next_bp;
					State ns = State(x,e);
					ns.free_edges++;
					if (ns.free_edges<=MAX_FREE && ns.score>10e-14) {
						pq.push(ns);
					}
				}
				fit--;
				if (fit!=bps.begin() && !forward) {
					fit--;
					pair<int, unsigned int> prev_bp=*fit;
					e.first=other_bp;
					e.second=prev_bp;
					State ps = State(x,e);
					ps.free_edges++;
					if (ps.free_edges<=MAX_FREE && ps.score>10e-14) {
						pq.push(ps);
					}
				}
			}


			//extend using normal edges
			State next = x.extend_normal();
			if (next.score>0) {
				if (next.free_edges<=MAX_FREE && next.score>10e-14) {
					pq.push(next);
				}
			}

		}
		#pragma omp critical
		{
			//cerr << thread_id << ":" << i << endl;
			if (best_state.score>10e-14) {
				best_state.print();
			}
			if (best_cycle.score>10e-14) {
				cout << "*" << best_cycle.is_cycle(true);
				best_cycle.print();
			}
		}

	}

	return 0;
}
