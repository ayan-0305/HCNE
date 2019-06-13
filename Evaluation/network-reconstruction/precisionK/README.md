# GraphEvaluation

precK : cpp file to evaluate the all pair reconstruction of graph from generated embeddings
Run : g++ -std=c++11 precK.cpp -o b.out
	  ./b.out <embeddings file> <graph file> <threshold outputfile <precision outputfile> <jaccard outputfile> <sim techniue ("cos"/"euc")>

#formats of files
--------------------
embeddings file : 
{#nodes #dims
id1 vec1
id2 vec2
...
}
--------------------
graph file
{id1 id2
id2 id3
..
}
--------------------
precision outputfile
{f1 f2
..
}
returns a ranked list of node pairs with decreasing similarity where:
f1 is an integer denoting whether corresponding node pair is connected (1) or not (0)
f2 is a floating point number denoting the corresponding value of precision@K
--------------------
Sample graph file: karate.txt
Sample embedding file: karate-emb.txt
Sample precision@K output file: prec_karate.txt
