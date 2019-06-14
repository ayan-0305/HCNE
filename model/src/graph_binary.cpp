// File: graph_binary.cpp
// -- graph handling source
//-----------------------------------------------------------------------------
// Community detection
// Based on the article "Fast unfolding of community hierarchies in large networks"
// Copyright (C) 2008 V. Blondel, J.-L. Guillaume, R. Lambiotte, E. Lefebvre
//
// And based on the article
// Copyright (C) 2013 R. Campigotto, P. Conde Céspedes, J.-L. Guillaume
//
// This file is part of Louvain algorithm.
// 
// Louvain algorithm is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// Louvain algorithm is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with Louvain algorithm.  If not, see <http://www.gnu.org/licenses/>.
//-----------------------------------------------------------------------------
// Author   : E. Lefebvre, adapted by J.-L. Guillaume and R. Campigotto
// Email    : jean-loup.guillaume@lip6.fr
// Location : Paris, France
// Time	    : July 2013
//-----------------------------------------------------------------------------
// see README.txt for more details


#include <fstream>
#include "graph_binary.h"
#define NLINKS 10000000



Graph::Graph() {
  nb_nodes = 0;
  nb_links = 0ULL;
  
  total_weight = 0.0L;
  sum_nodes_w = 0;
}


//compute the maximum of three unsigned long
inline unsigned long max3(unsigned long a,unsigned long b,unsigned long c){
  a = (a > b) ? a : b;
  return (a > c) ? a : c;
}

Graph::Graph(char *filename, char *filename_w, int type) {
  unsigned long e1 = NLINKS;
  ifstream finput;
	
  nb_nodes = 0;
  nb_links = 0ULL;

  finput.open(filename, fstream::in);
  if (finput.is_open() != true) {
    cerr << "The file " << filename << " does not exist" << endl;
    exit(EXIT_FAILURE);
  }

  // read file as edge list
  vector<pair<unsigned long, unsigned long> > edges;
  edges.resize(e1);

  while (finput >> edges[nb_links].first >> edges[nb_links].second) {
    nb_nodes = max3(nb_nodes, edges[nb_links].first, edges[nb_links].second);
    nb_links++;
    if (nb_links == e1) {
      e1 += NLINKS;
      edges.resize(e1);
    }
  }

  nb_nodes++;
  edges.resize(nb_links);
  /*
  for (int i = 0; i < edges.size(); i++) {
    cerr << edges[i].first << "-"  << edges[i].second << " " ;
  }
  cerr << endl;
  */

  // convert edge list to sparse graph format
  unsigned long i,u,v;
  vector<unsigned long> d(nb_nodes, 0);
  
  for (unsigned int i = 0; i < nb_links; i++) {
    d[edges[i].first]++;
    d[edges[i].second]++;
  }
  
  /*  for (int i = 0; i < nb_nodes; i++) {
    cerr << d[i] << " " ;
  }
  cerr << endl;
  */
  degrees.resize(nb_nodes);
  degrees[0] = d[0];
  for (int i = 1; i < nb_nodes; i++) {
    degrees[i] = degrees[i - 1] + d[i];
    d[i - 1] = 0;
  }
  d[nb_nodes-1]=0;

  /*
  cerr << nb_links << endl;
  for (int i = 0; i < nb_nodes; i++) {
    cerr << degrees[i] << " ";
  }
  cerr << endl;
  */
  links.resize(2*nb_links);

  for (i = 0; i < nb_links; i++) {
    u = edges[i].first;
    v = edges[i].second;

    //    if (v==198)
    //if (degrees[u-1] + d[u] > 2990000 || degrees[v-1] + d[v] > 2990000)
    //    cerr << u << " " << degrees[u-1] << " " <<  d[u] << " " << v << " " << degrees[v-1] << " " << d[v] << endl;

    if (u>0) {
      links[degrees[u-1] + d[u]] = v;
    } else {
      links[d[u]] = v;
    }
    d[u]++;
    if (v>0) {
      links[degrees[v-1] + d[v]] = u;
    } else {
      links[d[v]] = u;
    }
    d[v]++;
  }

  // Compute total weight
  total_weight = 0;
  for (int i=0 ; i<nb_nodes ; i++)
    total_weight += (long double)weighted_degree(i);

  nb_links *= 2;
  nodes_w.assign(nb_nodes, 1);
  sum_nodes_w = nb_nodes;

  //  cerr << degrees[nb_nodes-1] << " " << nb_links << endl;
}

/*
Graph::Graph(char *filename, char *filename_w, int type) {
  ifstream finput;
  finput.open(filename,fstream::in | fstream::binary);
  if (finput.is_open() != true) {
    cerr << "The file " << filename << " does not exist" << endl;
    exit(EXIT_FAILURE);
  }

  // Read number of nodes on 4 bytes
  finput.read((char *)&nb_nodes, sizeof(int));
  if (finput.rdstate() != ios::goodbit) {
    cerr << "The file " << filename << " is not a valid graph" << endl;
    exit(EXIT_FAILURE);
  }
  
  // Read cumulative degree sequence: 8 bytes for each node
  // cum_degree[0]=degree(0); cum_degree[1]=degree(0)+degree(1), etc.
  degrees.resize(nb_nodes);
  finput.read((char *)&degrees[0], nb_nodes*sizeof(unsigned long long));

  // Read links: 4 bytes for each link (each link is counted twice)
  nb_links = degrees[nb_nodes-1];
  links.resize(nb_links);
  finput.read((char *)(&links[0]), nb_links*sizeof(int));

  // IF WEIGHTED, read weights: 10 bytes for each link (each link is counted twice)
  weights.resize(0);
  total_weight = 0.0L;
  if (type==WEIGHTED) {
    ifstream finput_w;
    finput_w.open(filename_w,fstream::in | fstream::binary);
    if (finput_w.is_open() != true) {
      cerr << "The file " << filename_w << " does not exist" << filename << endl;
      exit(EXIT_FAILURE);
    }

    weights.resize(nb_links);
    finput_w.read((char *)(&weights[0]), nb_links*sizeof(long double));
    if (finput_w.rdstate() != ios::goodbit) {
      cerr << "The file " << filename_w << " does not correspond to valid weights for the graph" << filename << endl;
      exit(EXIT_FAILURE);
    }
  }

  // Compute total weight
  for (int i=0 ; i<nb_nodes ; i++)
    total_weight += (long double)weighted_degree(i);
  
  nodes_w.assign(nb_nodes, 1);
  sum_nodes_w = nb_nodes;
}
*/

Graph
Graph::subgraph(vector<int> nodes) {
  Graph sub;

  sub.nb_nodes = nodes.size();
  sub.degrees.resize(nodes.size());
  sub.weights.resize(0);
  sub.total_weight = 0.0L;

  vector<unsigned int> innodes(nb_nodes, 0);
  for (unsigned int i = 0; i < nodes.size(); i++) {
    innodes[nodes[i]] = i + 1;
  }

  for (unsigned int node = 0; node<nodes.size(); node++) {
    sub.degrees[node] = (node > 0)?sub.degrees[node - 1]:0;

    pair<vector<int>::iterator, vector<long double>::iterator > p = neighbors(nodes[node]);

    for (int i=0 ; i<nb_neighbors(nodes[node]) ; i++) {
      if (innodes[*(p.first + i)] != 0) {
	//if (weights.size()!=0)
	//  cout << " (" << *(p.first+i) << " " << *(p.second+i) << ")";
	//else
	
	sub.degrees[node]++;
	sub.nb_links++;
	sub.links.push_back(innodes[*(p.first + i)]-1);

	if (weights.size() != 0) {
	  sub.weights.push_back(*(p.second + i));
	  sub.total_weight += *(p.second + i);
	}
      }
    }
  }

  // Compute total weight
  for (int i=0 ; i<sub.nb_nodes ; i++)
    sub.total_weight += (long double)sub.weighted_degree(i);

  sub.nodes_w.assign(sub.nb_nodes, 1);
  sub.sum_nodes_w = sub.nb_nodes;

  return sub;
}

long double
Graph::max_weight() {
  long double max = 1.0L;

  if (weights.size()!=0)
    max = *max_element(weights.begin(),weights.end());
  
  return max;
}

void
Graph::assign_weight(int node, int weight) {
  sum_nodes_w -= nodes_w[node];

  nodes_w[node] = weight;

  sum_nodes_w += weight;
}

void
Graph::add_selfloops() {
  vector<unsigned long long> aux_deg;
  vector<int> aux_links;

  unsigned long long sum_d = 0ULL;

  for (int u=0 ; u < nb_nodes ; u++) {
    pair<vector<int>::iterator, vector<long double>::iterator> p = neighbors(u);
    int deg = nb_neighbors(u);

    for (int i=0 ; i < deg ; i++) {
      int neigh = *(p.first+i);
      aux_links.push_back(neigh);
    }

    sum_d += (unsigned long long)deg;

    if (nb_selfloops(u) == 0.0L) {
      aux_links.push_back(u); // add a selfloop
      sum_d += 1ULL;
    }

    aux_deg.push_back(sum_d); // add the (new) degree of vertex u
  }

  links = aux_links;
  degrees = aux_deg;
  
  nb_links += (unsigned long long)nb_nodes;
}

void
Graph::display() {
  for (int node=0 ; node<nb_nodes ; node++) {
    pair<vector<int>::iterator, vector<long double>::iterator > p = neighbors(node);
    cout << node << ":" ;
    for (int i=0 ; i<nb_neighbors(node) ; i++) {
      if (weights.size()!=0)
	cout << " (" << *(p.first+i) << " " << *(p.second+i) << ")";
      else
	cout << " " << *(p.first+i);
    }
    cout << endl;
  }
}

void
Graph::display_reverse() {
  for (int node=0 ; node<nb_nodes ; node++) {
    pair<vector<int>::iterator, vector<long double>::iterator > p = neighbors(node);
    for (int i=0 ; i<nb_neighbors(node) ; i++) {
      if (node>*(p.first+i)) {
	if (weights.size()!=0)
	  cout << *(p.first+i) << " " << node << " " << *(p.second+i) << endl;
	else
	  cout << *(p.first+i) << " " << node << endl;
      }
    }
  }
}

bool
Graph::check_symmetry() {
  int error = 0;
  for (int node=0 ; node<nb_nodes ; node++) {
    pair<vector<int>::iterator, vector<long double>::iterator > p = neighbors(node);
    for (int i=0 ; i<nb_neighbors(node) ; i++) {
      int neigh = *(p.first+i);
      long double weight = *(p.second+i);

      pair<vector<int>::iterator, vector<long double>::iterator > p_neigh = neighbors(neigh);
      for (int j=0 ; j<nb_neighbors(neigh) ; j++) {
	int neigh_neigh = *(p_neigh.first+j);
	long double neigh_weight = *(p_neigh.second+j);

	if (node==neigh_neigh && weight!=neigh_weight) {
	  cout << node << " " << neigh << " " << weight << " " << neigh_weight << endl;
	  if (error++==10)
	    exit(0);
	}
      }
    }
  }
  return (error==0);
}

void
Graph::display_binary(char *outfile) {
  ofstream foutput;
  foutput.open(outfile ,fstream::out | fstream::binary);

  foutput.write((char *)(&nb_nodes),sizeof(int));
  foutput.write((char *)(&degrees[0]),sizeof(unsigned long long)*nb_nodes);
  foutput.write((char *)(&links[0]),sizeof(int)*nb_links);
}
