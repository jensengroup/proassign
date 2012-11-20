from time import time

def lines(filename):				#Returns list_[lines][elements] from the source data.
	data = open(filename, 'r')
	#list_ = tuple()
	list_=[]
	for line in data:
		list_.append(line.rsplit())
		#list_+=(line,)
	data.close()
	return list_

list_=lines('1CEX_shiftx2.txt')

def histogram_element(list_):
	final=[]
	count = 0
	for element in range(len(list_[0])):
		hist={}
		for line in list_:
			s=line[element]
			if s in hist:
				hist[s] += 1
			else:
				hist[s] = 1
		v=list(hist.values())
		k=list(hist.keys())
		final.append(int(k[v.index(max(v))]))
	for i in range(len(final)):
		if final[i]==i:
			count+=1
	return final, count, len(final)

def histogram_structure(list_):
	hist={}
	for line in list_:
		s=line
		if s in hist:
			hist[s] += 1
		else:
			hist[s] = 1
	v=list(hist.values())
	k=list(hist.keys())
	final= k[v.index(max(v))]
	return final

A=time()
x,y,z= histogram_element(list_)


print x
print y,z,float(y)/z
print time()-A
#x=histogram_structure(list_)
#print x
