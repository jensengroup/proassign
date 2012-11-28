from time import time
from numpy import genfromtxt, arange, exp
import random
from copy import deepcopy
import pylab

def lines(filename):				#Returns list_[lines][elements] from the source data.
	data = open(filename, 'r')
	list_ = []
	for line in data:
		list_.append(line.rsplit())
	data.close()
	return list_

def str_(filename):			#outputs a temporary file of {(res.no, CS_type): CS}, with LESS unnecesary data :O
	list_=lines(filename)
	d={}
	structure=[]
	for line in range(len(list_)):
		if line==0:
			d[(line,list_[line][7])]= float(list_[line][10])
			structure.append(list_[line][6])
			l=0
		else:
			if list_[line][5] == list_[line-1][5]:
				d[(l,list_[line][7])]= float(list_[line][10])
			else:
				l+=1
				d[(l,list_[line][7])]= float(list_[line][10])
				structure.append(list_[line][6])
		if list_[line][6]=='GLY' and  list_[line][7]=='HA3':#Have no idea what the HA entry in cs files is, since there is HA, HA2 AND HA3 values listed. So I use HA3 as HA for GLY.
			d[l,'HA']= d[l,'HA3']
			d[l,'CB']='N/A'
		if list_[line][6]=='PRO':
			d[l,'H']='N/A'
			d[l,'N']='N/A'
	return d, structure, l+1

def histogram_(name):
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
	return final


def histogram(B):
	count = 0
	for i in range(len(B)):
		if B[i]==i:
			count+=1
	return count, len(B), float(count)/len(B)
	
def getdE(i,j,data, B, y):	
	dE=-(data[i][B[j]]+data[j][B[i]]-data[i][B[i]] - data[j][B[j]]) / y
	return dE



A=['1CEX']#['1A2P','1CEX','1DMB','1GZI','1HCB','1LDE', '1LLW','1RLC','1UBQ','2GB1','2YP2','2ZKR','3S9I','3TTQ','3VUO','4AMC','4D94','4HOL']
s=1000000

for l in A:
	name=l+'_camshift_real_NA.txt'
	data=genfromtxt(name)
	x=len(data[0])
	y=sum(data[0])	
	B=arange(x)
	random.shuffle(B)
	B_min=deepcopy(B)
	E_tot=0
	c_acc=0
	c_rej=0
	E_=[]
	E_min=0
	for k in range(s):
		i, j = random.randint(0,x-1), random.randint(0,x-1)
		dE=getdE(i,j,data, B, y)
		if dE<=0:
			B[i], B[j] = B[j], B[i]
			E_tot+=dE
		else:
			p=random.random()
			if p<exp(-dE*10000): #6
				B[i], B[j] = B[j], B[i]
				E_tot+=dE
				c_acc+=1
			else:
				c_rej+=1
		if E_tot < E_min:
			E_min=E_tot
			B_min=deepcopy(B)
		E_.append(E_tot)
	#print c_acc, c_rej
	#print l, histogram(B_min)#, B_min
	'''pylab.plot(E_, marker='.', linestyle='None')
	#pylab.ylim(min(E_)-abs(min(E_))*1.5,abs(max(E_[100*d_le:]))*1.5)
	pylab.savefig('hist_test')'''
	hist_=histogram_(name)



a,b,c = str_('1CEX.str')
structure={}
for i in range(len(b)):
	structure[i]=b[i]

output=open('output_camshift.txt', 'w')

for i in range(len(b)):
	A=str(i)+'\t'+structure[i]+'\t %6.2f'% (100*float(data[i][i])/y)+'%'
	B=str(B_min[i])+'\t'+structure[B_min[i]]+'\t %6.2f'% (100*float(data[i][B_min[i]])/y)+'%'
	C=str(hist_[i])+'\t'+structure[hist_[i]]+'\t %6.2f'% (100*float(data[i][hist_[i]])/y)+'%'
	output.write(A+'\t'*3+B+'\t'*3+C+'\n')
