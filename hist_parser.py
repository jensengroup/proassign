from time import time
from numpy import genfromtxt

def histogram(name):
	final=[]
	count = 0
	data=genfromtxt(name)
	S=sum(data[0])
	for i in range(len(data[0])):
		hist={}
		for j in range(len(data[0])):
				hist[j]=data[i][j] #have double checked that it is not [j][i].
		v=list(hist.values())
		k=list(hist.keys())
		final.append(int(k[v.index(max(v))]))
	for i in range(len(final)):
		if final[i]==i:
			count+=1
	return final, count, len(final)

A=[ '1A2P', '1CEX', '1DMB','1GZI', '1HCB', '1ubq', '2gb1']
for i in A:
	x,y,z= histogram(i+'_shiftx2.txt')
	print i, y,z,float(y)/z
	for l in range(z):
		if l not in x:
			1#print 'Error in most probable assignment: Multiples of same value. \n Value not found: %d' %l
