from skmultilearn.cluster.networkx import NetworkXLabelGraphClusterer
from skmultilearn.ensemble import LabelSpacePartitioningClassifier
from skmultilearn.cluster import LabelCooccurrenceGraphBuilder
import networkx as nx
import re
import numpy as np
import pickle
import sys
from sklearn.ensemble import RandomForestClassifier
from skmultilearn.problem_transform import BinaryRelevance,ClassifierChain,LabelPowerset
from sklearn.naive_bayes import GaussianNB
from sklearn.metrics.pairwise import cosine_similarity
from sklearn import metrics
from sklearn.model_selection import cross_val_predict
from sklearn.linear_model import LogisticRegressionCV, LogisticRegression
from scipy.spatial import distance
from sklearn.cross_validation import KFold, cross_val_score
from sklearn.cross_validation import StratifiedKFold
from sklearn.cross_validation import train_test_split
from sklearn.metrics import f1_score, accuracy_score, precision_recall_fscore_support, adjusted_rand_score
from sklearn.metrics.cluster import normalized_mutual_info_score, homogeneity_score 
from sklearn import tree
from sklearn.cluster import KMeans

def purity_score(y_true, y_pred):
    y_true=np.array(y_true)
    y_pred=np.array(y_pred)
    y_voted_labels = np.zeros(y_true.shape)
    labels = np.unique(y_true)
    ordered_labels = np.arange(labels.shape[0])
    for k in range(labels.shape[0]):
        y_true[y_true==labels[k]] = ordered_labels[k]
    labels = np.unique(y_true)
    bins = np.concatenate((labels, [np.max(labels)+1]), axis=0)
    for cluster in np.unique(y_pred):
        hist, _ = np.histogram(y_true[y_pred==cluster], bins=bins)
        winner = np.argmax(hist)
        y_voted_labels[y_pred==cluster] = winner
    return accuracy_score(y_true, y_voted_labels)

nlabels=sys.argv[1]
d=pickle.load(open(nlabels,"rb"))
l=[]
for k in d.keys():
    for x in d[k]:
        l.append(int(x))

l=list(set(l))
l=sorted(l)
print max(l),min(l)
d1={}
for k in d.keys():
    d1[k]=[]
    for i in range(0,len(l)):
        d1[k].append(0)
    for j in range(0,len(d[k])):
        x=d[k][j]
        y=l.index(x)
        d1[k][y]=1

D={}
embedding=sys.argv[2]
f=open(embedding,'r')
c=0
for line in f:
    c+=1
    if c==1:
        continue
    t=line.split()
    if len(t)==0:
        continue
    D[t[0]]=[]
    for i in range(1,len(t)):
        D[t[0]].append(float(t[i]))

key=d.keys()

for i in range(0,100):
    np.random.shuffle(key)

acc=[]
nmi=[]
ari=[]
hom=[]
wt=[]
for i in range(0,len(l)):
    L=[]
    L1=[]

    for k in key:
        try:
            L.append(D[k])
        except:
            continue
        L1.append(d1[k][i])

    X=np.array(L)
    Y=np.array(L1)

    kmeans = KMeans(n_clusters=len(l), random_state=0).fit(X)
    #print id1
    g=kmeans.labels_.tolist()



    fa=accuracy_score(g,L1)
    fb=normalized_mutual_info_score(g,L1)
    fc=adjusted_rand_score(g,L1)
    fd=homogeneity_score(g,L1)
    fe=f1_score(g,L1,average='weighted')
    acc.append(fa)
    nmi.append(fb)
    ari.append(fc)
    hom.append(fd)
    wt.append(fe)
    print fa, fb

print "Accuracy=",np.mean(acc),"NMI=",np.mean(nmi)
print "Variances=",np.var(acc),np.var(nmi)