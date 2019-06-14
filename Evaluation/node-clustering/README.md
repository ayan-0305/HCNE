For evaluation of generated embedding on node clustering, run the command:

python node-cluster.py node-label embedding-file

where node-label -> file is a .p pickle file containing a python dictionary whose key is node id and value is a list of node labels embedding-file -> denotes the file containing embedding vectors for each node in the format

   {id1 vec1 
   id2 vec2 
   ... }

The program outputs the NMI scores for node clustering. 
