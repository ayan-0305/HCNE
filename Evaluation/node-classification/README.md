For multi-label classification using generated embedding vectors, run the command:

python multi-label.py node-label embedding-file

where  node-label -> file is a .p pickle file containing a python dictionary whose key is node id and value is a list of node labels
       embedding-file -> denotes the file containing embedding vectors for each node in the format
       
       {id1 vec1 
       id2 vec2 
       ... }
       
The program outputs the Micro-F1 and Macro-F1 scores for multi-label node classification.
