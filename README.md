# LouvainNE
LouvainNE: Hierarchical Louvain Method for High Quality and Scalable Graph Embedding

The 'model' folder contains the codebase for generating the HCNE embedding vectors for nodes in an input graph.
To generate the embedding vectors, run the following command inside the 'model' folder to generate embedding vectors for nodes in the karate club network:
./louvain karate.txt -e 0.005 -a 20.0 -n 20 > karate-emb.txt

where karate.txt -> edgelist file for the karate network
      -e 0.005 sets the precision limit of Louvain
      -a 20.0 is the 1/alpha parameter for combining the embeddings
      -n 20 sets the dimension of the embedding vectors to 20
      karate-emb.txt -> output embedding file containing embedding vectors of nodes in the format:
      {id1 vec1
       id2 vec2
       ...
      }


The 'Evaluation' folder contains the code for evaluating generated embedding vectors on different downstream prediction tasks as well as generating a sampled graph representative of the original graph. 

For our experiments, we are using the following datasets:

Blogcatalog: http://socialcomputing.asu.edu/datasets/BlogCatalog3
youtube: http://socialnetworks.mpi-sws.org/data-imc2007.html
flickr: http://socialnetworks.mpi-sws.org/data-imc2007.html
