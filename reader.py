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
		try:
			for i in range(2,len(list_[2])):
				if float(list_[line][i])!=0:
					d.setdefault((int(list_[line][0])-1,list_[0][i]), list_[line][i])
				else:
					d.setdefault((int(list_[line][0])-1,list_[0][i]), 'N/A')
		except IndexError:
			pass
	return d


def str_(filename):			#outputs a temporary file of {(res.no, CS_type): CS}, with LESS unnecesary data :O
	list_=lines(filename)
	d={}
	structure=[]
	for line in range(len(list_)):
		if line==0:
			d.setdefault((line,list_[line][7]), float(list_[line][10]))
			structure.append(list_[line][6])
			l=0
		else:
			if list_[line][5] == list_[line-1][5]:
				d.setdefault((l,list_[line][7]), float(list_[line][10]))
			else:
				l+=1
				d.setdefault((l,list_[line][7]), float(list_[line][10]))
				structure.append(list_[line][6])
	return d, structure, l+1


def cs_(filename):
	list_=lines(filename)
	d={}
	for line in range(len(list_)):
		if line==0:
			l=0
			d.setdefault((l,list_[line][3]), float(list_[line][4]))
		else:
			if list_[line][1] == list_[line-1][1]:
				d.setdefault((l,list_[line][3]), float(list_[line][4]))
			else:
				l+=1
				d.setdefault((l,list_[line][3]), float(list_[line][4]))
	return d
