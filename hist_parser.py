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
				hist[j]=data[j][i] #Dunno if (i,j) or (j,i) if it matters at all.
		v=list(hist.values())
		k=list(hist.keys())
		final.append(int(k[v.index(max(v))]))
	for i in range(len(final)):
		if final[i]==i:
			count+=1
	return final, count, len(final)

x,y,z= histogram('1CEX_shiftx2.txt')


print x
print y,z,float(y)/z

