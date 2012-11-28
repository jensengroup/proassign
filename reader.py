from numpy import *

### Most of the N/A statements are in here in addition to the pair function in assign.py to be able to force gly / pro into place by e.g. putting '-100' instead. But non-complete CS pairs in the bmrb database messes that idea up.

def lines(filename):				#Returns list_[lines][elements] from the source data.
	data = open(filename, 'r')
	list_ = []
	for line in data:
		list_.append(line.rsplit())
	data.close()
	return list_


def cam_(filename):			#outputs a temporary file of {(res.no, CS_type): CS}
	list_=lines(filename)
	d={}
	for line in range(2,len(list_)):
		for i in range(2,len(list_[2])):
			l=list_[0][i]
			if l=='HN':
				l='H'
			if float(list_[line][i])!=0:
				d[(line-2,l)]= float(list_[line][i])
			else:
				d[(line-2,l)]= 'N/A'
		if list_[line][1]=='GLY' and float(list_[line][2])!=0: #second criteria prevents that the first and last element in camshift gets edited if it's all zero.
			d[(line-2,'CB')]= 'N/A'		#sets CB to 'N/A'.
		if list_[line][1]=='PRO' and float(list_[line][2])!=0: #float due to camshift being unconsistant with these things for some reason.
			d[(line-2,'H')]= 'N/A'
			d[(line-2,'N')]= 'N/A'
	return d, line-1


def str_(filename):			#outputs a temporary file of {(res.no, CS_type): CS}.
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
		if list_[line][6]=='GLY' and  list_[line][7]=='HA3':# Forces HA3 as HA for GLY in the BMRB file.
			d[l,'HA']= d[l,'HA3']
			d[l,'CB']='N/A'
		if list_[line][6]=='PRO':
			d[l,'H']='N/A'
			d[l,'N']='N/A'
	return d, structure, l+1


def cs_(filename):
	list_=lines(filename)
	d={}
	for line in range(len(list_)):
		if line==0:
			l=0
			d[l,list_[line][3]]= float(list_[line][4])
		else:
			if list_[line][1] == list_[line-1][1]:
				d[l,list_[line][3]]= float(list_[line][4])
			else:
				l+=1
				d[l,list_[line][3]]= float(list_[line][4])
		if list_[line][2]=='GLY' and  list_[line][3]=='HA3':#Have no idea what the HA entry in cs files is, since there is HA, HA2 AND HA3 values listed. So I use HA3 as HA for GLY.
			d[l,'HA']= d[l,'HA3']
			d[l,'CB']='N/A'
		if list_[line][2]=='PRO':
			d[l,'H']='N/A'
			d[l,'N']='N/A'
	return d, l+1
