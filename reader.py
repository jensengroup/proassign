from numpy import *

def lines(filename):				#Returns list_[lines][elements] from the source data.
	data = open(filename, 'r')
	list_ = []
	for line in data:
		list_.append(line.rsplit())
	data.close()
	return list_


def cam_(filename):			#outputs a temporary file of {(res.no, CS_type): CS}, with some unnecesary data.
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
			d[(line-2,'CB')]= -100		#sets CB to -100. This makes sure that all GLY gets paired together.
		if list_[line][1]=='PRO' and float(list_[line][2])!=0: #float due to camshift being unconsistant with these things for some reason.
			print list_[line][2]
			d[(line-2,'H')]= -100
			d[(line-2,'N')]= -100
	return d, line-1


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
			d[l,'CB']=-100
		if list_[line][6]=='PRO':
			d[l,'H']=-100
			d[l,'N']=-100
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
			d[l,'CB']=-100
		if list_[line][2]=='PRO':
			d[l,'H']=-100
			d[l,'N']=-100
	return d, l+1
