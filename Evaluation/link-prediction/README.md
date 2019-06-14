Steps for performing link prediction:

1. For obtaining the pruned graph by removing a fraction of edges from original graph for link prediction,
run the command: python prunegraph.py original-graph-file pruned-graph-file

The edgelist for the pruned graph will be saved in pruned-graph-file

2. Generate embeddings for all nodes using the pruned-graph-file as the input graph.

3. Perform link prediction to report the accuracy on the generated embeddings in Step 2 using the command:
python link-pred.py original-graph-file pruned-graph-file pruned-graph-embedding-file 
where pruned-graph-embedding-file contains the embeddings learned on the pruned graph.

