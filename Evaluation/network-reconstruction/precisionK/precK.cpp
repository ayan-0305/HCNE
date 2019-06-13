#include <bits/stdc++.h>
#include <unistd.h>

#define vi vector<int>
#define vvi vector<vi>
#define vd vector<double>
#define vvd vector<vd>
#define pii pair<int,int> 
#define pb push_back
#define ff first
#define ss second

using namespace std;

int print_bool = 1;
map<int, int> id_map, inv_id_map;

int get_intersection_card(vi &first, vi &second){
	vi intersect(first.size());
	auto it = set_intersection (first.begin(), first.end(), second.begin(), second.end(), intersect.begin());
	return it - intersect.begin();
}

//assumes first and second received are of equal sizes
//returns jaccard, normalised
pair<double, double> Set_sim(vi &first, vi &second){
	sort(first.begin(), first.end());
	sort(second.begin(), second.end());
	int int_card = get_intersection_card(first, second);
	double vaules0 =  (double)int_card / (2.0*first.size() - int_card);
	double values1 = (double)int_card / first.size();

	return make_pair(vaules0, values1);
} 

double cosine_sim(vd &first, vd &second){
	int n = first.size();
	double dot_prod=0, mod1=0, mod2=0;
	for(int i=0;i<n;i++){
		double a,b;
		a = first[i]; b = second[i];
		dot_prod +=  a*b;
		mod1 += a*a;
		mod2 += b*b; 
	}
	return dot_prod/(sqrt(mod1)*sqrt(mod2));
}

double euclid_sim(vd &first, vd &second){
	int n = first.size();
	double sqr_dist=0;
	for(int i=0;i<n;i++){
		sqr_dist += ((first[i]-second[i])*(first[i]-second[i]));
	}
	return -1.0*sqrt(sqr_dist);
}

double get_sim(vd &first, vd &second, string sim_crit = "cos"){
	if(sim_crit=="cos")
		return cosine_sim(first, second);
	else if(sim_crit=="euc")
		return euclid_sim(first, second);
	return 0.0;
}
void check_missing(string filename){
	ifstream myfile;
	myfile.open(filename);
	vi done(id_map.size(), 0);
	int node_count, dim;
	myfile >> node_count;
	if(node_count != id_map.size()){
		cout << node_count << " " << id_map.size() << endl;
		cout << "DANGER!!1-=-=-=-=-=-=-=-=-=-=-=-=-=-=-DANGER!!1-=-=-=-=-=-=-=-=-=-=-=-=-=-=-DANGER!!1-=-=-=-=-=-=-=-=-=-=-=-=-=-=-DANGER!!1-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n\n\n\n\n\n";
	}
	myfile >> dim;
	cout << dim << endl;
	cout << "checking missing node" << endl;
	cout << "Reading embeddings from input file..." << filename << "with " << node_count << " embeddings" << endl;
	for(int i=0 ; i<node_count ; i++){
		if((i%1000==0 && print_bool))
			cout << i << endl;
		int node_id;
		myfile >> node_id;
		done[id_map[node_id]]=1;
		double sum = 0;
		for(int j=0;j<dim;j++){
				double temp;
				myfile >> temp;
		}
	}

	for(int i=0;i<id_map.size();i++){
		if(done[i]==0)
			cout << "BC : "<< i << " " << inv_id_map[i] << endl;
 	}
	cout << endl;
	myfile.close();
}
void load_embeddings(vvd &emb, string filename){
	ifstream myfile;
	myfile.open(filename);
	int node_count, dim;
	myfile >> node_count;
	if(node_count != id_map.size()){
		cout << node_count << " " << id_map.size() << endl;
		cout << "DANGER!!1-=-=-=-=-=-=-=-=-=-=-=-=-=-=-DANGER!!1-=-=-=-=-=-=-=-=-=-=-=-=-=-=-DANGER!!1-=-=-=-=-=-=-=-=-=-=-=-=-=-=-DANGER!!1-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n\n\n\n\n\n";
	}
	myfile >> dim;
	cout << dim << endl;
	cout << "Allocating memory...\n";
	emb.resize(node_count, vd(dim));
	cout << "Reading embeddings from input file..." << filename << "with " << node_count << " embeddings" << endl;
	for(int i=0 ; i<node_count ; i++){
		if((i%1000==0 && print_bool))
			cout << i << endl;
		int node_id;
		myfile >> node_id;
		double sum = 0;
		for(int j=0;j<dim;j++){
				myfile >> emb[id_map[node_id]][j];
		}
	}
	cout << endl;
	myfile.close();
}

int load_graph(set<pii> &graph, vvi &adj_list, string filename){
	
	cout << "Creating edgelist from input file...\n";
	int count = 0 ;
	ifstream myfile;
	myfile.open(filename);
	while(myfile){
		count++;
		if(count%10000==0 && print_bool)
			cout << count << endl;
		
		int a,b;
		myfile >> a >> b;
		if(!myfile){
			continue;
		}
		if(id_map.find(a) == id_map.end()){
			int k = id_map.size();
			id_map[a] = k;
			inv_id_map[k] = a;
		}
		if(id_map.find(b) == id_map.end()){
			int k = id_map.size();
			id_map[b] = k;
			inv_id_map[k] = b;
		}
		a = id_map[a];
		b = id_map[b];
		if(a > b)
			swap(a,b);
		graph.insert(make_pair(a,b));
	}
	cout << endl;
	myfile.close();

	cout << "Creating adjlist from loaded edgelist... " << id_map.size() << endl;
	int max_node = id_map.size()-1;	
	adj_list.resize(max_node+1,vi());
	cout << "yes\n\n\n";
	count = 0 ;
	for(auto it = graph.begin(); it!=graph.end(); it++){
		count++;
		if(count%1000==0 && print_bool)
			cout << count << endl;
		int a,b;
		a = it->first;
		b = it->second;
		adj_list[a].pb(b);
		adj_list[b].pb(a);
	}
	cout << "preprocessing adjlist..." << endl;
	for(int nd = 0; nd<=max_node; nd++){
		if(nd%1000==0 && print_bool)
			cout << nd << endl;
		sort(adj_list[nd].begin(), adj_list[nd].end());
	}
	return count;
}

double get_threshold_at_rank(vvd &emb, int rank, string sim_crit="cos"){
	cout << "Getting Threshold at rank " << rank << endl;
	vd sims;
	int n = emb.size();
	for(int i=0;i<n;i++){
		for(int j=i+1;j<n;j++){
			sims.push_back(get_sim(emb[i], emb[j], sim_crit));
		}
	}

	sort(sims.begin(), sims.end());

	return sims[sims.size() - rank];
}

//redundant
double get_max_sims(vvd &emb, string sim_crit){
	double possible_max = (double)emb[0].size();

	double maxi = DBL_MAX;
	int n = emb.size();
	for(int i=0;i<n;i++){
		for(int j=i+1;j<n;j++){
			maxi = max(maxi, get_sim(emb[i], emb[j], sim_crit));
		}
	}
	return maxi;
}

// 0:tp, 1:fn, 2:fp, 3:tn 
void evaluate(vvd &emb, set<pii> &graph, double threshold, vi &accuracy, string sim_crit){
	int n = emb.size();
	int dim = emb[0].size();
	int count = 0;
	for(int i=0;i<n;i++){
		for(int j=i+1;j<n;j++){
			int index = 0;
			if(graph.find(make_pair(i,j)) == graph.end()){
				index |= 2;
			}
			if(get_sim(emb[i], emb[j], sim_crit) < threshold){
				index |= 1;
			}
			accuracy[index]++;

			count++;
			if(count%10000000 == 0 && print_bool)
				cout << count << endl;
		}
	}
}

void print_precision_at_k(vvd &emb, set<pii> &graph, string filename, string sim_crit="cos"){
	clock_t t0, t1;
	t0 = clock();
	cout << "gRaph size = " << graph.size() << endl;
	cout << "Evaluating the precision @ k..." << endl;
	vector<pair<double, pair<int,int> > > sim_id;
	int n = emb.size();
	for(int i=0;i<n;i++){
		if(i%1000==0)
			cout << i << " " << sim_id.size()<< endl;
		for(int j=i+1;j<n;j++){
			sim_id.push_back(make_pair(-1.0 * get_sim(emb[i], emb[j], sim_crit), make_pair(i,j)));
		}
	}
	cout << "sorting the similarities..." << endl;
	sort(sim_id.begin(), sim_id.end());
	
	int crt, wrg;
	crt = wrg = 0;
	ofstream myfile;
	myfile.open(filename);
	cout << "writing the precision @ k..." << endl;
	for(int i=0;i<sim_id.size();i++){
		if(i%10000000==0)
			cout << i << endl;
		if(graph.find(make_pair(sim_id[i].ss.ff,sim_id[i].ss.ss)) == graph.end()){
			wrg++;
			myfile << 0 << " ";
		}
		else{
			crt++;
			myfile << 1 << " ";
		}
		myfile << (double)crt/(crt+wrg) << endl;
		if(crt > 0.5*graph.size())
			break;
	}
	myfile.close();
	t1 = clock();
	cout << "\n Done. Time Taken : " << ((double)(t1 - t0) / CLOCKS_PER_SEC) << endl;
	cout << endl;
}
void print_set_based_accuracy(vvi &adj_list, vvd &emb, string filename, string sim_crit="cos", int k=-1)
{
	clock_t t0, t1;
	t0 = clock();

	cout << "Evaluating the Jaccard Coefficients..." << endl;
	ofstream myfile;
	myfile.open(filename);
	myfile << "Node\t\t\t\tDegree\tJaccardCoef\tNormalisedCoef\n";
	cout << "Node\tDegree\tJaccardCoef\tNormalisedCoef\n";
	int n = adj_list.size();
	cout << n << " nodes in the graph" << endl;
	for(int src=0; src<n; src++){
		if(src%1000==0 && print_bool)
			cout << src << endl;
		vector<pair<double, int> > sim_id(n, make_pair(0.0,0));
		int degree = adj_list[src].size();
		for(int dst=0;dst<n;dst++){
			if(dst==src)
				continue;
			sim_id[dst] = make_pair(-1.0*get_sim(emb[src], emb[dst], sim_crit), dst);
		}
		sort(sim_id.begin(), sim_id.end());
		vi pred_nbr(degree);
		for(int i=0;i<degree;i++){
			pred_nbr[i] = sim_id[i].second;
		}
		auto sim_set = Set_sim(pred_nbr, adj_list[src]);
		myfile << setw(6) << inv_id_map[src] << "\t"<< setw(10) << degree << "\t"<< setw(10) << sim_set.first << "\t"<< setw(10) << sim_set.second << endl;
	}
	myfile.close();
	t1 = clock();
	cout << "\n Done. Time Taken : " << ((double)(t1 - t0) / CLOCKS_PER_SEC) << endl;
	cout << endl;
}
void print_accuracy_for_threshold(double threshold, vvd &emb, set<pii> &graph, string filename, string sim_crit, int flag_first=0){
	clock_t t0, t1;
	t0 = clock();

	cout << "Evaluating for the threshold " << threshold << endl;

	// 0:tp, 1:fn, 2:fp, 3:tn 
	vi accuracy(4,0);
	evaluate(emb, graph, threshold, accuracy, sim_crit);
	
	ofstream myfile;
	if(flag_first){
		myfile.open(filename);
		myfile << "threshold\t:\ttp\t\tfn\t\tfp\t\t\ttn\t\tAccuracy\tPrecision\tRecall\t\tTime" << endl;
	}
	else{
		myfile.open(filename, ios_base::app);
	}
	if(sim_crit=="euc")
		threshold = -1.0*threshold;
	cout << setw(8) << threshold << "\t:\t";
	myfile << threshold << "\t:\t";
	for(int i=0;i<accuracy.size();i++){
		cout << accuracy[i] << "\t";
		myfile << setw(8) << accuracy[i] << "\t";
	}
	//accuracy
	cout << ((double)accuracy[0] + accuracy[3])/(accuracy[0] + accuracy[1] + accuracy[2] + accuracy[3]) << "\t";
	myfile << setw(12) << ((double)accuracy[0] + accuracy[3])/(accuracy[0] + accuracy[1] + accuracy[2] + accuracy[3]) << "\t";
	//precision
	cout << ((double)accuracy[0])/(accuracy[0] + accuracy[2]) << "\t";
	myfile << setw(12) << ((double)accuracy[0])/(accuracy[0] + accuracy[2]) << "\t";
	//recall
	cout << ((double)accuracy[0])/(accuracy[0] + accuracy[1]) << "\t";
	myfile << setw(12) << ((double)accuracy[0])/(accuracy[0] + accuracy[1]) << "\t";
	
	t1 = clock();
	myfile << setw(12) << ((double)(t1 - t0) / CLOCKS_PER_SEC) << "s" << endl;
	myfile.close();
	cout << "\n Done. Time Taken : " << ((double)(t1 - t0) / CLOCKS_PER_SEC) << endl;
	cout << endl;
}
int main(int argc, char** argv){

	// string emb_filename = "../emb/pgp_512.emd";
	// string graph_filename = "../graph/pgp.edgelist";
	// string outfilename = "./output/AccOutput_512.txt";
	string emb_filename(argv[1]);
	string graph_filename(argv[2]);
	string outfilename(argv[3]);
	string precKfilename(argv[4]);
	string Jacfilename(argv[5]);
	string sim_measure(argv[6]);
	cout << argv[1] << " " << argv[2] << " " << argv[3] << " " << argv[4] << " " << argv[5] << endl; 

	set<pii> graph;
	vvi adj_list;
	int edge_count = load_graph(graph, adj_list, graph_filename);
	
	// check_missing(emb_filename);
	vvd emb;
	load_embeddings(emb, emb_filename);
	
	double min_sim;
	if(sim_measure == "euc")
		min_sim = -2.0 * sqrt(emb[0].size());

	//print_set_based_accuracy(adj_list, emb, Jacfilename, sim_measure);
	print_precision_at_k(emb, graph, precKfilename, sim_measure);
	for(double th = -1.0; th <=1; th+=0.1){
		double req_th = th;
		if(sim_measure=="euc"){
			req_th = ((-th+1.0)/2)*min_sim;
			cout << min_sim << ", required threshold ->" << req_th << endl;
		}
		print_accuracy_for_threshold(req_th, emb, graph, outfilename, sim_measure, th<=-0.95);
	}

	// print_accuracy_for_threshold(get_threshold_at_rank(emb, edge_count), emb, graph, outfilename);

	return 0;
}