import networkx as nx
import re
import sys
import numpy as np
import pickle
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
from sklearn.metrics import f1_score, accuracy_score, precision_recall_fscore_support
from sklearn import tree
from sklearn.utils import shuffle


nlabels=sys.argv[1]
d=pickle.load(open(nlabels,"rb"))
l=[]
for k in d.keys():
    print d[k]
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
thres=.3
c=0
for line in f:
    c+=1
    t=line.split()
    if len(t)<=2:
        continue
    D[t[0]]=[]
    for i in range(1,len(t)):
        D[t[0]].append(float(t[i]))

key=d.keys()

for i in range(0,100):
    np.random.shuffle(key)

mic=0
mac=0
while(True):  
    L=[]
    L1=[]

    for k in d1.keys():
        try:
            L.append(D[k])
        except:
            continue
        L1.append(d1[k])

    X=np.array(L)
    Y=np.array(L1)
    X_s,Y_s=shuffle(X,Y)
    size=[0.2]
    Mic=[]
    Mac=[]
    Wt=[]
    Acc=[]
    if mic>=thres:
        break
    for j in range(0,len(size)):
        X_train, X_test, Y_train, Y_test = train_test_split(X_s, Y_s, test_size=size[j])
        #k_fold = KFold(len(Y), n_folds=10, shuffle=True, random_state=0)
        clf = ClassifierChain(LogisticRegression())
        #clf = tree.DecisionTreeClassifier()
        #clf=RandomForestClassifier()
        clf.fit(X_train,Y_train)
        Y_predicted = clf.predict(X_test)
        Mic.append(f1_score(Y_test, Y_predicted, average='micro'))
        Mac.append(f1_score(Y_test, Y_predicted, average='macro'))
        Wt.append(f1_score(Y_test, Y_predicted, average='weighted'))
        Acc.append((accuracy_score(Y_test, Y_predicted)))
        mic+=Mic[j]
        mac+=Mac[j]
print "Micro-F1=",float(mic)/float(len(size))*1.0
print "Macro-F1=",float(mac)/float(len(size))*1.0
