import numpy as np
from sklearn.metrics import roc_auc_score
from scipy.spatial import distance
import time
from sklearn import metrics

f=open("karate.txt",'r')
edges=[]
nodes=[]
pose=[]

for line in f:
	t=line.split()
	pose.append((t[0],t[1]))
	edges.append((t[1],t[0]))
	nodes.append(t[0])
	nodes.append(t[1])

nodes=list(set(nodes))

k=len(pose)

itera=0
noe=[]
for i in range(0,len(nodes)):
	for j in range(1,len(nodes)):
		noe.append((nodes[i],nodes[j]))
		print len(noe)
		if len(noe)>=3*len(pose):
			print "breaking Here"
			#time.sleep(3)
			break
	if len(noe)>=3*len(pose):
		break

noe=list(set(noe).difference(set(pose)))
noe=list(set(noe).difference(set(edges)))


edges1=edges[:k]
pose1=pose[:k]
print len(noe),len(pose1)
AUC=[]

while (itera<=100):

	itera+=1
	print "iter=",itera
	f1=open("karate-emb.txt",'r')

	for i in range(0,100):
		np.random.shuffle(nodes)
		np.random.shuffle(edges)
		np.random.shuffle(pose)
		np.random.shuffle(noe)


	#np.random.shuffle(pose)

	d={}
	c=0
	for line in f1:
		t=line.split()
		c+=1
		if len(t)==2:
			continue
		d[t[0]]=[]
		for i in range(1,len(t)):
			d[t[0]].append(float(t[i]))


	print len(d)

	L=[]
	L1=[]
	L2=[]
	for i in range(0,len(pose1)):
		m=pose1[i][0]
		m1=pose1[i][1]
		l=noe[i][0]
		l1=noe[i][1]
		#print d1[m],d1[m1],d1[l],d1[l1]
		try:
			s=0.0
			s1=0.0
			for j in range(0,len(d[m])):
				s+=d[m][j]*d[m1][j]
				s1+=d[l][j]*d[l1][j]
			#print s,s1
			#s=distance.euclidean(d[m],d[m1])
			#s1=distance.euclidean(d[l],d[l1])
			L2.append((2,s))
			L2.append((1,s1))
		except:
			continue
	sorted(L2, key=lambda x: (x[1],x[1]), reverse=True)

	for i in range(0,len(L2)):
		L.append(L2[i][0])
		L1.append(L2[i][1])

	y = np.array(L)
	pred = np.array(L1)
	fpr, tpr, thresholds = metrics.roc_curve(y, pred, pos_label=2)
	print "auc=",metrics.auc(fpr,tpr)
	AUC.append(metrics.auc(fpr,tpr))

print np.mean(AUC)
