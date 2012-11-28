import re

cam_name='1A2P.cam'
str_name='1A2P.str'
cs_name='1A2P.cs'

[ '1A2P', '1CEX', '1DMB','1GZI', '1HCB', '1ubq', '2gb1']

def cam_(filename):	
	data = open(filename, 'r')
	list_ = []
	for line in data:
		m_obj= re.search(r"[A-Z]{3}", line)
		if m_obj:
			list_.append(m_obj.group(0))
	return list_

def str_(filename):
	data = open(filename, 'r')
	list_ = []
	lal=[]
	for line in data:
		list_.append(line.rsplit())
	a=0
	for i in list_:
		if i[5]==a:
			pass
		else:
			lal.append(i[6])
			a=i[5]
	return lal

def cs_(filename):
	data = open(filename, 'r')
	list_ = []
	lal=[]
	for line in data:
		list_.append(line.rsplit())
	a=0
	for i in list_:
		if i[1]==a:
			pass
		else:
			lal.append(i[2])
			a=i[1]
	return lal




def compare():
	C=cam_(cam_name)
	A=str_(str_name)
	B=cs_(cs_name)
	print 'str','cs','cam'
	for i in range(len(A)):
		print A[i], B[i], C[i], A[i]==B[i], A[i]==C[i], i
	print len(A), len(B), len(C)

compare()
