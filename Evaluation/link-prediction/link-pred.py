import networkx as nx
import re
import sys
import numpy as np
import pickle
from sklearn.metrics.pairwise import cosine_similarity
from sklearn import metrics
from sklearn.model_selection import cross_val_predict
from sklearn.linear_model import LogisticRegressionCV, LogisticRegression
from scipy.spatial import distance

def hadamard(l1,l2):
    l3=[]
    for i in range(0,len(l1)):
        l3.append(float(l1[i]*l2[i]))
    return l3

def average(l1,l2):
    l3=[]
    for i in range(0,len(l1)):
        l3.append((l1[i]+l2[i])*0.5)
    return l3

def wt_L1(l1,l2):
    l3=[]
    for i in range(0,len(l1)):
        l3.append(abs(l1[i]-l2[i]))
    return l3

def wt_L2(l1,l2):
    l3=[]
    for i in range(0,len(l1)):
        l3.append(abs(l1[i]-l2[i])*abs(l1[i]-l2[i]))
    return l3

graph=sys.argv[1]
prunedgraph=sys.argv[2]
f=open(prunedgraph,'r')
f1=open(graph,'r')

pose=[]
poset=[]
nege=[]
noe=[]
nege1=[]
poe=[]

for line in f:
    t=line.split()
    pose.append((t[0],t[1]))
    pose.append((t[1],t[0]))
    nege1.append((t[0],t[1]))
    

for line in f1:
	t=line.split()
	poset.append((t[0],t[1]))
	nege.append((t[0],t[1]))
	poset.append((t[1],t[0]))


pose=list(set(pose))

embedding=sys.argv[3]

d={}

f=open(embedding,'r')
c=0
for line in f:
    c+=1
    t=line.split()
    if len(t)==2:
        continue
    d[t[0]]=[]
    for i in range(1,len(t)):
        d[t[0]].append(float(t[i]))

key=d.keys()

np.random.shuffle(pose)
np.random.shuffle(key)


for i in range(0,len(key)):
    for j in range(1,len(key)):
        noe.append((key[i],key[j]))
        #noe=list(set(noe).difference(set(poset)))
        print len(noe)
        if len(noe)>=3*len(poset):
            break
    if len(noe)>=3*len(poset):
        break

noe=list(set(noe).difference(set(poset)))
poe=list(set(nege).difference(set(nege1)))

#nege1=list(set(nege).difference(set(pose)))


np.random.shuffle(noe)

t=len(poe)

noe=noe[:t]


L=[]
L1=[]


for i in range(0,len(poe)):
    k=poe[i][0]
    k1=poe[i][1]
    h=average(d[k],d[k1])

    L.append(h)
    L1.append(1)
    #print len(L)

    k=noe[i][0]
    k1=noe[i][1]

    h=average(d[k],d[k1])
    L.append(h)
    L1.append(0)
    print len(L)


X=np.array(L)
Y=np.array(L1)
logreg=LogisticRegressionCV()
predicted = cross_val_predict(logreg, X, Y, cv=10)
print metrics.accuracy_score(Y, predicted)
print metrics.classification_report(Y, predicted)
print metrics.roc_auc_score(Y, predicted)
