import sys
graph=sys.argv[1]
prunedgraph=sys.argv[2]

frac_removed=0.4

f=open(graph,'r')
import numpy as np
d={}
for line in f:
	t=line.split()
	d[t[0]]=0
	d[t[1]]=0

f=open(graph,'r')
for line in f:
	t=line.split()
	d[t[0]]+=1
	d[t[1]]+=1

c=0
f=open(graph,'r')
f1=open(prunedgraph,'w')
for line in f:
	t=line.split()
	if np.random.uniform(0,1)>=(1-frac_removed) and d[t[0]]>1 and d[t[1]]>1:
		d[t[0]]-=1
		d[t[1]]-=1
		continue
	f1.write(str(t[0]))
	f1.write(" ")
	f1.write(str(t[1]))
	f1.write("\n")

