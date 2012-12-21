import random
from numpy import exp, sqrt, pi, log, average, std, zeros, savetxt, mean
import pylab
import copy
from time import time, clock
import create_spins_from_bmrb

def lines(filename):				#Returns list_[lines][elements] from the source data.
	data = open(filename, 'r')
	list_ = []
	for line in data:
		list_.append(line.rsplit())
	data.close()
	return list_


def read(filename):
	list_=lines(filename)
	d=[]
	for line in range(2,len(list_)):
		l=[]
		for i in range(1,len(list_[line])):
			l.append(float(list_[line][i]))
		d.append(l)
	return d

def sequence_(filename):
	sequence=[]
	list_=lines(filename)
	for i in list_:
		for l in i:
			for n in l:
				sequence.append(n)
	return sequence


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
			d[(line-2,'HA3')]=d[(line-2,'HA3')]
			del d[(line-2,'HA')]
		if list_[line][1]=='PRO' and float(list_[line][2])!=0: #float due to camshift being unconsistant with these things for some reason.
			d[(line-2,'H')]= 'N/A'
			d[(line-2,'N')]= 'N/A'
	return d,line-1


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
	return d,l+1


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
		if list_[line][2]=='GLY' and  list_[line][3]=='HA':#Have no idea what the HA entry in cs files is, since there is HA, HA2 AND HA3 values listed. 
			d[l,'HA']= 'N/A'
	return d,l+1



def AA_(filename):			#outputs a temporary file of {(res.no, CS_type): CS}.
	list_=lines(filename)
	d={}
	structure=[]
	for line in range(len(list_)):
		if line==0:
			d[(line,list_[line][5])]= float(list_[line][8])
			structure.append(list_[line][4])
			l=0
		else:
			if list_[line][3] == list_[line-1][3]:
				d[(l,list_[line][5])]= float(list_[line][8])
			else:
				l+=1
				d[(l,list_[line][5])]= float(list_[line][8])
				structure.append(list_[line][4])
	return d,l+1



def HSQC_make(HSQC):
	t=[]
	for i in HSQC:
		t_=[]
		for l in range(len(CS_pairs)):
			t_.append('N/A')
		for l in range(len(CS_pairs)):
			if CS_pairs[l]=='N':
				t_.append(i[0])
			elif CS_pairs[l]=='H':
				t_.append(i[1])
			else:
				t_.append('N/A')
		t.append(t_)
	return t

def HNCO_make(HNCO):
	t=[]
	for i in HNCO:
		t_=[]
		for l in range(len(CS_pairs)):
			if CS_pairs[l]=='C':
				t_.append(i[0])
			else:
				t_.append('N/A')
		for l in range(len(CS_pairs)):
			if CS_pairs[l]=='N':
				t_.append(i[1])
			elif CS_pairs[l]=='H':
				t_.append(i[2])
			else:
				t_.append('N/A')
		t.append(t_)
	return t

def HNcoCA_make(HNcoCA):
	t=[]
	for i in HNcoCA:
		t_=[]
		for l in range(len(CS_pairs)):
			if CS_pairs[l]=='CA':
				t_.append(i[0])
			else:
				t_.append('N/A')
		for l in range(len(CS_pairs)):
			if CS_pairs[l]=='N':
				t_.append(i[1])
			elif CS_pairs[l]=='H':
				t_.append(i[2])
			else:
				t_.append('N/A')
		t.append(t_)
	return t

def HNcoHA_make(HNcoHA):
	t=[]
	for i in HNcoHA:
		t_=[]
		for l in range(len(CS_pairs)):
			if CS_pairs[l]=='HA':
				t_.append(i[0])
			else:
				t_.append('N/A')
		for l in range(len(CS_pairs)):
			if CS_pairs[l]=='N':
				t_.append(i[1])
			elif CS_pairs[l]=='H':
				t_.append(i[2])
			else:
				t_.append('N/A')
		t.append(t_)
	return t


steps=300001 #maxsteps to be run

Names=['2L15'] 
theory_input='BMRB' # 'BMRB', 'shiftx2' or 'camshift'
CS_pairs=['H', 'HA', 'HA2', 'HA3', 'N', 'C', 'CA', 'CB']
CS_sigma_cam={'N':1.12, 'CA':0.44, 'CB':0.52, 'C':0.98, 'HA':0.12,'HA2':0.12,'HA3':0.12, 'H':0.17}#{'H':0.56,'HA':0.28,'HA2':0.28,'HA3':0.28,'N':3.01,'C':1.38, 'CA':1.3, 'CB':1.36}		#dictionary with value=variance for key=CS type. #Kohlhoff 2009
CS_sigma_cs={'N':1.12, 'CA':0.44, 'CB':0.52, 'C':0.98, 'HA':0.12,'HA2':0.12,'HA3':0.12, 'H':0.17} 	#shiftx2, CD and CG variance set to same as similar shift types, Han, B  Liu, Y, 2011
CS_sigma_xp={'N':0.1, 'CA':0.1, 'CB':0.1, 'C':0.1, 'HA':0.02,'HA2':0.02,'HA3':0.02, 'H':0.02}
CS_sigma_AA={'N':0.1, 'CA':0.1, 'CB':0.1, 'C':0.1, 'HA':0.02,'HA2':0.02,'HA3':0.02, 'H':0.02}

spin_systems=['HSQC','HNCO','HNcoCA']#,'HNcoHA'] #HA for gly
sequence=sequence_('cspa/sequence.txt')

move_single_percentage=1



#######################################################################

def spin_make():
	t=[]
	for i in spin_systems:
		if i=='HSQC':
			try:
				HSQC=read('cspa/HSQC.list')
				t.append(HSQC_make(HSQC))
			except NameError:
				pass
		if i=='HNCA':
			try:
				HNCA=read('cspa/HNCA.list')
				t.append(HNCA_make(HNCA))
			except NameError:
				pass
		if i == 'HNCACB':
			try:
				HNCACB=read('cspa/HNCACB.list')
				t.append(HNCACB_make(HNCACB))
			except NameError:
				pass
		if i == 'HNCO':
			try:
				HNCO=read('cspa/HNCO.list')
				t.append(HNCO_make(HNCO))
			except NameError:
				pass
		if i == 'HNcoCA':
			try:
				HNcoCA=read('cspa/HNcoCA.list')
				t.append(HNcoCA_make(HNcoCA))
			except NameError:
				pass
		if i == 'HNcoCACB':
			try:
				HNcoCACB=read('cspa/HNcoCACB.list')
				t.append(HNcoCACB_make(HNcoCACB))
			except NameError:
				pass
		if i=='HNcoHA':
			try:
				HNcoHA=read('cspa/HNcoHA.list')
				t.append(HNcoHA_make(HNcoHA))
			except NameError:
				pass
		if i =='HNHA':
			try:
				HNHA=read('cspa/HNHA.list')
				t.append(HNHA_make(HNHA))
			except NameError:
				pass
	return t

#HSQC=[N(i),H(i)]
#HNCA=[CA(i(+),i-1(-)),N(i), H(i)]
#HNCACB=[CAB(i,i-1), N(i), H(i)]
#HNCO=[C(i-1), N(i), H(i)]
#HNcoCA=[CA(i-1), N(i), H(i)]
#HNcaCACB=[CAB(i-1),N(i),H(i)]
#HNcoHA=?[HA(i-1),N(i),H(i)]
#HNHA=[HA(i(+),i-1(-)),N(i),H(i)]


if theory_input=='camshift':		# Determines which sigma dictionary to use
	CS_sigma=CS_sigma_cam
elif theory_input=='shiftx2':
	CS_sigma=CS_sigma_cs
elif theory_input=='BMRB':
	CS_sigma=CS_sigma_xp
elif theory_input=='AA':
	CS_sigma=CS_sigma_AA



def get_average(t,i,T):
	t_=[]
	for m in range(len(CS_pairs)):
		temp=[]
		for l in t:
			if i!=0:
				try:
					temp.append(float(l[i-1][len(CS_pairs)+m]))
				except ValueError:	
					pass
			try:
				temp.append(float(l[i][m]))
			except ValueError:	
				pass
		if len(temp)!=0:
			t_.append(sum(temp)/float(len(temp)))
		else:
			t_.append('N/A')
	for m in range(len(CS_pairs)):
		temp=[]
		for l in t:
			try:
				temp.append(float(l[i][len(CS_pairs)+m]))
			except ValueError:	
				pass
			if i!=len(l)-1:
				try:
					temp.append(float(l[i+1][m]))
				except ValueError:	
					pass
		if len(temp)!=0:
			t_.append(sum(temp)/float(len(temp)))
		else:
			t_.append('N/A')	
	T.append(t_)
	return T


	

def getE(t,d,r): #calculates total energy difference of the two systems. 
	E=0
	T=[]
	for m in r:
		t_=[]
		for l in range(len(t)):
			for n in range(len(t)):
				if n<l:
					E+=getE_part(t[l][m],t[n][m],CS_sigma_xp)
				if m!=0:
					E+=getE_part(t[l][m-1][len(CS_pairs):],t[n][m][:-len(CS_pairs)],CS_sigma_xp)
				if m!=len(d)-1:
					E+=getE_part(t[l][m][len(CS_pairs):],t[n][m+1][:-len(CS_pairs)],CS_sigma_xp)
		T=get_average(t,m,T)
	for i in range(len(r)):
		s=d_make(d,r[i])
		E+=getE_part(T[i],s)
	return E

def d_make(d,i):
	if i==0:
		s=d[i][:]
		for l in range(len(d[i])):
			s=['N/A']+s
	else:
		s=d[i-1][:]+d[i][:]
	return s




def pair(t,d,d_le):		#t = a set of CS types to be used
	d_=[]	
	for i in range(d_le):
		temp=[]
		for l in t:
			try:
				temp.append(d[(i, l)])
			except KeyError:
				temp.append('N/A')
		d_.append(temp)
	return d_

def max_length(t,d):
	l=0
	for i in t:
		if len(i)>l:
			l=len(i)
	if len(d)>l:
		l=len(d)
	return l

def set_length(t,d,le):
	T=[]
	for l in t:
		t_=l
		if len(l)<le:
			count=le-len(l)	
			for n in range(count):
				temp=[]
				for i in range(len(l[0])):
					temp.append('N/A')
				t_.append(temp)
		T.append(t_)
	if len(d)<le:
		count=le-len(d)
		for n in range(count):
			temp=[]
			for i in range(len(d[0])):
				temp.append('N/A')
			d.append(temp)
	return T,d

def getdE_part(t,d,i,j,CS_sigma=CS_sigma):
	dE=getE_part(t[i],d[j],CS_sigma)+getE_part(t[j],d[i],CS_sigma)-getE_part(t[i],d[i],CS_sigma)-getE_part(t[j],d[j],CS_sigma)
	return dE

def getE_part(t,d,CS_sigma=CS_sigma):		#get the current energy difference of t and d_e at site i.
	dE=0
	CS=CS_pairs*2
	#print t
	for n in range(len(t)):
		try:
			dc=t[n]-d[n]
			dE += (dc**2) /(2 * CS_sigma[CS[n]]**2)
		except TypeError:	#pass if there's one or both of the probed chemical shifts that's missing.
			pass
	return dE

def determine_bias_spin(t):
	count=0
	for n in t:
		if type(n) == float:
			count=1
		else:
			pass
	if count==1:
		bias=True
	else:
		bias=False
	return bias

def determine_bias(t):
	bias=[]
	for n in t:
		if type(n) == float:
			bias.append(True)
		else:
			bias.append(False)
	return bias

'''def getdE_change_HA(t,d,i,j,k):
	p=random.randint(2,3)
	dE=0
	if spin_systems[k]=='HNHA' or spin_systems[k]=='HNcoHA':
		for n in range(len(CS_pairs)):
			if CS_pairs[n]=='HA':
				for l in range(len(CS_pairs)):
					if CS_pairs[l]='HA2':
						for m in range(len(t)):
							for m!=HNcoHA_index:
								if t[m][i][n]!='N/A' and t[HNcoHA_index][i][n]!='N/A' #eliminating the bias will make dE=0 for every HA->HA2 move.
								dE+=getE_part(t[m][i],t[HNcoHA_index][i],CS_sigma_xp)
						
						r[k][i]=t[k][i][:]
						r[k][i][l], r[k][i][j]=r[k][i][l], 'N/A'
						getE_part(r[k][i]'''



def getdE_move(t,d,i,j,k):		#get energy difference for interchanging i and j
	dE=0
	try:
		if sequence[i]!='G':
			if sequence[j]!='G':
				dE=getdE_move_single(t,d,i,j,k)
				return dE
			else:
				i,j=j,i
		else:
			if sequence[j]!='G':
				pass
			else:
				1
	except IndexError:
		pass
	
	return 

def getdE_move_single(t,d,i,j,k):		#get energy difference for interchanging i and j
	if i==j:
		return 0, 0		#no move.
	dE=0
	dB=0
	bias=[]
	T=[]
	s=[]
	d_le=len(d)-1
	for l in range(len(spin_systems)):	#noting if any float values exist in the l[i],l[j],l[i-1],l[j-1],l[i+1].l[j+1] spinsystems
		bias.append(determine_bias_spin(t[l][i]))
		bias.append(determine_bias_spin(t[l][j]))
		if i!=0:
			bias.append(determine_bias_spin(t[l][i-1]))
		else:
			bias.append(False)
		if j!=0:
			bias.append(determine_bias_spin(t[l][j-1]))
		else:
			bias.append(False)
		if i!=d_le:
			bias.append(determine_bias_spin(t[l][i+1]))
		else:
			bias.append(False)
		if j!=d_le:
			bias.append(determine_bias_spin(t[l][i-1]))
		else:
			bias.append(False)

	for l in range(len(spin_systems)):
		if bias[6*k]==bias[6*k+1]:	#No bias if both or neither of l[i],l[j] exist.
			if bias[6*k]==False:
				pass	#No energy if neither of above exist.
			else:
				if l!=k:
					dE+=getdE_part(t[l],t[k],i,j,CS_sigma_xp)
				# the below prevents indexerrors.
				if i!=0:
					dE-=getE_part(t[l][i-1][len(CS_pairs):],t[k][i][:-len(CS_pairs)],CS_sigma_xp)
					dE+=getE_part(t[l][i-1][len(CS_pairs):],t[k][j][:-len(CS_pairs)],CS_sigma_xp)
				if j!=0:
					dE-=getE_part(t[l][j-1][len(CS_pairs):],t[k][j][:-len(CS_pairs)],CS_sigma_xp)
					dE+=getE_part(t[l][j-1][len(CS_pairs):],t[k][i][:-len(CS_pairs)],CS_sigma_xp)

				if i!=d_le:
					dE-=getE_part(t[l][i+1][:-len(CS_pairs)],t[k][i][len(CS_pairs):],CS_sigma_xp)
					dE+=getE_part(t[l][i+1][:-len(CS_pairs)],t[k][j][len(CS_pairs):],CS_sigma_xp)
				if j!=d_le:
					dE-=getE_part(t[l][j+1][:-len(CS_pairs)],t[k][j][len(CS_pairs):],CS_sigma_xp)
					dE+=getE_part(t[l][j+1][:-len(CS_pairs)],t[k][i][len(CS_pairs):],CS_sigma_xp)		

		else:
			if l!=k:
				dE+=getdE_part(t[l],t[k],i,j,CS_sigma_xp)
				if bias[6*l]!=bias[6*l+1]:	#The number of energyterms added and substracted should remain constant, which the below enforces.
					if bias[6*k]==True and bias[6*l]==True:
						dB+=getE_part(t[l][i],t[k][i],CS_sigma_xp)
					if bias[6*k]==True and bias[6*l]==False:
						dB-=getE_part(t[l][j],t[k][i],CS_sigma_xp)
					if bias[6*k]==True and bias[6*l+1]==True:
						dB+=getE_part(t[l][j],t[k][j],CS_sigma_xp)
					if bias[6*k]==True and bias[6*l+1]==False:
						dB-=getE_part(t[l][i],t[k][j],CS_sigma_xp)

			# the below prevents indexerrors.
			if i!=0:
				dE-=getE_part(t[l][i-1][len(CS_pairs):],t[k][i][:-len(CS_pairs)],CS_sigma_xp)
				dE+=getE_part(t[l][i-1][len(CS_pairs):],t[k][j][:-len(CS_pairs)],CS_sigma_xp)
			if j!=0:
				dE-=getE_part(t[l][j-1][len(CS_pairs):],t[k][j][:-len(CS_pairs)],CS_sigma_xp)
				dE+=getE_part(t[l][j-1][len(CS_pairs):],t[k][i][:-len(CS_pairs)],CS_sigma_xp)

			if i!=d_le:
				dE-=getE_part(t[l][i+1][:-len(CS_pairs)],t[k][i][len(CS_pairs):],CS_sigma_xp)
				dE+=getE_part(t[l][i+1][:-len(CS_pairs)],t[k][j][len(CS_pairs):],CS_sigma_xp)
			if j!=d_le:
				dE-=getE_part(t[l][j+1][:-len(CS_pairs)],t[k][j][len(CS_pairs):],CS_sigma_xp)
				dE+=getE_part(t[l][j+1][:-len(CS_pairs)],t[k][i][len(CS_pairs):],CS_sigma_xp)	
			# many o bias was had.
			if bias[6*l+2]==True and bias[6*l+4]==True and bias[6*l+3]==True and bias[6*l+5]==False and i!=d_le:
				dB+=getE_part(t[k][i][len(CS_pairs):],t[l][i+1][:-len(CS_pairs)],CS_sigma_xp)
			if bias[6*l+2]==True and bias[6*l+4]==True and bias[6*l+3]==False and bias[6*l+5]==True and i!=0:
				dB+=getE_part(t[l][i-1][len(CS_pairs):],t[k][i][:-len(CS_pairs)],CS_sigma_xp)
			if bias[6*l+2]==True and bias[6*l+4]==False and bias[6*l+3]==True and bias[6*l+5]==True and j!=d_le:
				dB-=getE_part(t[k][i][len(CS_pairs):],t[l][j+1][:-len(CS_pairs)],CS_sigma_xp)
			if bias[6*l+2]==False and bias[6*l+4]==True and bias[6*l+3]==True and bias[6*l+5]==True and j!=0:
				dB-=getE_part(t[l][j-1][len(CS_pairs):],t[k][i][:-len(CS_pairs)],CS_sigma_xp)
			if bias[6*l+2]==False and bias[6*l+4]==False and bias[6*l+3]==True and bias[6*l+5]==True:
				if j!=0:
					dB-=getE_part(t[l][j-1][len(CS_pairs):],t[k][i][:-len(CS_pairs)],CS_sigma_xp)
				if j!=d_le:
					dB-=getE_part(t[k][i][len(CS_pairs):],t[l][j+1][:-len(CS_pairs)],CS_sigma_xp)
			if bias[6*l+2]==True and bias[6*l+4]==True and bias[6*l+3]==False and bias[6*l+5]==False:
				if i!=0:
					dB+=getE_part(t[l][i-1][len(CS_pairs):],t[k][i][:-len(CS_pairs)],CS_sigma_xp)
				if i!=d_le:
					dB+=getE_part(t[k][i][len(CS_pairs):],t[l][i+1][:-len(CS_pairs)],CS_sigma_xp)
			if bias[6*l+2]==True and bias[6*l+4]==False and bias[6*l+3]==False and bias[6*l+5]==False and i!=0:
				dB+=getE_part(t[l][i-1][len(CS_pairs):],t[k][i][:-len(CS_pairs)],CS_sigma_xp)
			if bias[6*l+2]==False and bias[6*l+4]==True and bias[6*l+3]==False and bias[6*l+5]==False and i!=d_le:
				dB+=getE_part(t[k][i][len(CS_pairs):],t[l][i+1][:-len(CS_pairs)],CS_sigma_xp)
			if bias[6*l+2]==False and bias[6*l+4]==False and bias[6*l+3]==True and bias[6*l+5]==False and j!=0:
				dB-=getE_part(t[l][j-1][len(CS_pairs):],t[k][i][:-len(CS_pairs)],CS_sigma_xp)
			if bias[6*l+2]==False and bias[6*l+4]==False and bias[6*l+3]==False and bias[6*l+5]==True and j!=d_le:
				dB-=getE_part(t[k][i][len(CS_pairs):],t[l][j+1][:-len(CS_pairs)],CS_sigma_xp)
##i->j
			if bias[6*l+3]==True and bias[6*l+5]==True and bias[6*l+2]==True and bias[6*l+4]==False and j!=d_le:
				dB+=getE_part(t[k][j][len(CS_pairs):],t[l][j+1][:-len(CS_pairs)],CS_sigma_xp)
			if bias[6*l+3]==True and bias[6*l+5]==True and bias[6*l+2]==False and bias[6*l+4]==True and j!=0:
				dB+=getE_part(t[l][j-1][len(CS_pairs):],t[k][j][:-len(CS_pairs)],CS_sigma_xp)
			if bias[6*l+3]==True and bias[6*l+5]==False and bias[6*l+2]==True and bias[6*l+4]==True and i!=d_le:
				dB-=getE_part(t[k][j][len(CS_pairs):],t[l][i+1][:-len(CS_pairs)],CS_sigma_xp)
			if bias[6*l+3]==False and bias[6*l+5]==True and bias[6*l+2]==True and bias[6*l+4]==True and i!=0:
				dB-=getE_part(t[l][i-1][len(CS_pairs):],t[k][j][:-len(CS_pairs)],CS_sigma_xp)
			if bias[6*l+3]==False and bias[6*l+5]==False and bias[6*l+2]==True and bias[6*l+4]==True:
				if i!=0:
					dB-=getE_part(t[l][i-1][len(CS_pairs):],t[k][j][:-len(CS_pairs)],CS_sigma_xp)
				if i!=d_le:
					dB-=getE_part(t[k][j][len(CS_pairs):],t[l][i+1][:-len(CS_pairs)],CS_sigma_xp)
			if bias[6*l+3]==True and bias[6*l+5]==True and bias[6*l+2]==False and bias[6*l+4]==False:
				if j!=0:
					dB+=getE_part(t[l][j-1][len(CS_pairs):],t[k][j][:-len(CS_pairs)],CS_sigma_xp)
				if j!=d_le:
					dB+=getE_part(t[k][j][len(CS_pairs):],t[l][j+1][:-len(CS_pairs)],CS_sigma_xp)
			if bias[6*l+3]==True and bias[6*l+5]==False and bias[6*l+2]==False and bias[6*l+4]==False and j!=0:
				dB+=getE_part(t[l][j-1][len(CS_pairs):],t[k][j][:-len(CS_pairs)],CS_sigma_xp)
			if bias[6*l+3]==False and bias[6*l+5]==True and bias[6*l+2]==False and bias[6*l+4]==False and j!=d_le:
				dB+=getE_part(t[k][j][len(CS_pairs):],t[l][j+1][:-len(CS_pairs)],CS_sigma_xp)
			if bias[6*l+3]==False and bias[6*l+5]==False and bias[6*l+2]==True and bias[6*l+4]==False and i!=0:
				dB-=getE_part(t[l][i-1][len(CS_pairs):],t[k][j][:-len(CS_pairs)],CS_sigma_xp)
			if bias[6*l+3]==False and bias[6*l+5]==False and bias[6*l+2]==False and bias[6*l+4]==True and i!=d_le:
				dB-=getE_part(t[k][j][len(CS_pairs):],t[l][i+1][:-len(CS_pairs)],CS_sigma_xp)



	if j==i+1 or j==i-1:

		dE-=getE_part(t[k][i][len(CS_pairs):],t[k][i][:-len(CS_pairs)],CS_sigma_xp)
		dE-=getE_part(t[k][j][len(CS_pairs):],t[k][j][:-len(CS_pairs)],CS_sigma_xp)
		dE+=2*getE_part(t[k][j][len(CS_pairs):],t[k][i][:-len(CS_pairs)],CS_sigma_xp)	
		if bias[6*k]==True and bias[6*k+1]==False:
			if bias[6*k+2]==False:
				dB+=getE_part(t[k][i][len(CS_pairs):],t[k][i][:-len(CS_pairs)],CS_sigma_xp)
		if bias[6*k]==False and bias[6*k+1]==True:
			if bias[6*k+5]==False:
				dB+=getE_part(t[k][j][len(CS_pairs):],t[k][j][:-len(CS_pairs)],CS_sigma_xp)

	'''if j==i-1:
		B_E=dE
		B_B=dB
		dE-=getE_part(t[k][i][len(CS_pairs):],t[k][i][:-len(CS_pairs)],CS_sigma_xp)
		dE-=getE_part(t[k][j][len(CS_pairs):],t[k][j][:-len(CS_pairs)],CS_sigma_xp)
		dE+=2*getE_part(t[k][i][len(CS_pairs):],t[k][j][:-len(CS_pairs)],CS_sigma_xp)	
		if bias[6*k+1]==True and bias[6*k]==False:
			if bias[6*k+3]==False:
				dB+=getE_part(t[k][j][len(CS_pairs):],t[k][j][:-len(CS_pairs)],CS_sigma_xp)
		if bias[6*k+1]==False and bias[6*k]==True:
			if bias[6*k+6]==False:
				dB+=getE_part(t[k][i][len(CS_pairs):],t[k][i][:-len(CS_pairs)],CS_sigma_xp)
		A_E=dE
		A_B=dB
		dijE=A_E-B_E
		dijB=A_B-B_B'''



	t[k][i], t[k][j] = t[k][j], t[k][i]
	T=get_average(t,i,T)#Tij
	T=get_average(t,j,T)#Tji
	t[k][i], t[k][j] = t[k][j], t[k][i]
	T=get_average(t,i,T)#Tii
	T=get_average(t,j,T)#Tjj
	s.append(d_make(d,i))
	s.append(d_make(d,j))
	bias=[]
	bias.append(determine_bias(T[0]))
	bias.append(determine_bias(T[1]))
	bias.append(determine_bias(T[2]))
	bias.append(determine_bias(T[3]))
	bias.append(determine_bias(s[0]))
	bias.append(determine_bias(s[1]))
	dE+=getE_part(T[0],s[0])+getE_part(T[1],s[1])-getE_part(T[2],s[0])-getE_part(T[3],s[1])
	
	for l in range(len(CS_pairs)*2):
		p=[]
		p.append(bias[2][l]*bias[4][l])
		p.append(bias[3][l]*bias[5][l])
		p.append(bias[0][l]*bias[4][l])
		p.append(bias[1][l]*bias[5][l])
		if p[0]==True and p[1]==False and p[2]==False and p[3]==False:
			dB+=getE_part([T[2][l]],[s[0][l]])
		if p[0]==False and p[1]==True and p[2]==False and p[3]==False:
			dB+=getE_part([T[3][l]],[s[1][l]])
		if p[0]==False and p[1]==False and p[2]==True and p[3]==False:
			dB-=getE_part([T[0][l]],[s[0][l]])
		if p[0]==False and p[1]==False and p[2]==False and p[3]==True:
			dB-=getE_part([T[1][l]],[s[1][l]])
		if p[0]==True and p[1]==True and p[2]==False and p[3]==False:
			dB+=getE_part([T[2][l]],[s[0][l]])
			dB+=getE_part([T[3][l]],[s[1][l]])
		if p[0]==False and p[1]==False and p[2]==True and p[3]==True:
			dB-=getE_part([T[0][l]],[s[0][l]])
			dB-=getE_part([T[1][l]],[s[1][l]])
		if p[0]==False and p[1]==True and p[2]==True and p[3]==True:
			dB-=getE_part([T[0][l]],[s[0][l]])
		if p[0]==True and p[1]==False and p[2]==True and p[3]==True:
			dB-=getE_part([T[1][l]],[s[1][l]])
		if p[0]==True and p[1]==True and p[2]==False and p[3]==True:
			dB+=getE_part([T[2][l]],[s[0][l]])
		if p[0]==True and p[1]==True and p[2]==True and p[3]==False:
			dB+=getE_part([T[3][l]],[s[1][l]])
	return dE+dB, dE

def getdE_move_block(t,d,i,j):
	T=[]
	s=[]
	dE=0
	dB=0
	for m in range(len(spin_systems)):
		t[m][i], t[m][j] = t[m][j], t[m][i]
	T=get_average(t,i,T) #Tij
	T=get_average(t,j,T)#Tj
	for m in range(len(spin_systems)):
		t[m][i], t[m][j] = t[m][j], t[m][i]
	T=get_average(t,i,T) #Tii
	T=get_average(t,j,T)#Tjj

	s.append(d_make(d,i))
	s.append(d_make(d,j))

	bias=[]
	bias.append(determine_bias(T[0]))
	bias.append(determine_bias(T[1]))
	bias.append(determine_bias(s[0]))
	bias.append(determine_bias(s[1]))

	dE+=getE_part(T[0],s[1])+getE_part(T[1],s[0])-getE_part(T[0],s[0])-getE_part(T[1],s[1])
	for l in range(len(T[0])):
		if T[0][l] == True and d[0][l]==True and T[1][l]==False and d[1][l]==False:
			dB+=getE_part([T[0][l]],[d[0][l]])
		elif T[0][l] == False and d[0][l]==False and T[1][l]==True and d[1][l]==True:
			dB+=getE_part([T[1][l]],[d[1][l]])
		elif T[0][l] == True and d[0][l]==False and T[1][l]==False and d[1][l]==True:
			dB-=getE_part([T[0][l]],[d[1][l]])
		elif T[0][l] == False and d[0][l]==True and T[1][l]==True and d[1][l]==False:
			dB+=getE_part([T[1][l]],[d[0][l]])		
			
	return dE+dB, dE


if __name__ == '__main__':
	def run():
		for n in Names:

			if theory_input=='BMRB':
				d,d_le = str_(n+'.str')	#BMRB filename
			elif theory_input=='camshift':
				d,d_le = cam_(n+'.cam')	#camshift filename
			elif theory_input=='shiftx2':
				d,d_le = cs_(n+'.cs') # name of shiftx2 file
			elif theory_input=='AA':
				d,d_le = AA_(n+'_AA.str') # name of AA file
			else:
				print 'Error in theory_input'
				break

			d=pair(CS_pairs,d,d_le)

			#create_spins_from_bmrb.HSQC(d,CS_pairs)
			t=spin_make()
			le=max_length(t,d)
			t,d=set_length(t,d,le)
			for i in range(len(t)):
				random.shuffle(t[i])
			E_=[]
			L=range(len(t[0]))
			E_tot=getE(t,d,L)		
			T=[]
			for i in range(le):
				T=get_average(t,i,T)
			B=getE([T],d,L)
			print 0, E_tot, B
			count_1=0
			count_2=0
			count_3=0
			count_4=0
			for s in range(1,steps):
				n=random.random()
				i, j = random.randint(0,le-1), random.randint(0,le-1)
				
				if n<move_single_percentage:
					count_1+=1
					k=random.randint(0,len(t)-1)
					#if i!=j+1 and i!=j-1:
					dE,dE_tot=getdE_move_single(t,d,i,j,k)
					#else:
					#	dE=1
					if dE<=0:
						count_2+=1
						#E_before=getE(t,d,[i,j])
						t[k][i], t[k][j] = t[k][j], t[k][i]
						#E_after=getE(t,d,[i,j])
						#E_tot+=E_after-E_before
						E_tot+=dE_tot
						#if round(dE_tot,3)!=round(E_after-E_before,3):
						#	print dE_tot, E_after-E_before
					else:
						p=random.random()
						if p<exp(-dE/20):
							count_2+=1
							#E_before=getE(t,d,[i,j])
							t[k][i], t[k][j] = t[k][j], t[k][i]
							#E_after=getE(t,d,[i,j])
							#E_tot+=E_after-E_before
							E_tot+=dE_tot
				else:
					count_3+=1
					dE,dE_tot=getdE_move_block(t,d,i,j)
					if dE<=0:
						count_4+=1
						E_before=getE(t,d,[i,j])
						for m in range(len(spin_systems)):
							t[m][i], t[m][j] = t[m][j], t[m][i]
						E_after=getE(t,d,[i,j])
						#E_tot+=E_after-E_before
						E_tot+=dE_tot
						if round(dE_tot,3)!=round(E_after-E_before,3):
							print dE_tot, E_after-E_before
					else:
						p=random.random()
						if p<exp(-dE/20):
							count_4+=1
							#E_before=getE(t,d,[i,j])
							for m in range(len(spin_systems)):
								t[m][i], t[m][j] = t[m][j], t[m][i]
							#E_after=getE(t,d,[i,j])
							#E_tot+=E_after-E_before
							E_tot+=dE_tot
				E_.append(E_tot)
				#if s>10:
				#	if E_[-1]>E_[-2]:
				#		print 'Error'
				if s%5000==0:

					T=[]
					#E_t=getE(t,d,L)
					for i in range(le):
						T=get_average(t,i,T)
					B=getE([T],d,L)
					try:
						print s, E_tot, B#, float(count_2)/count_1, float(count_4)/count_3
					except ZeroDivisionError:
						pass
					count_1=0
					count_2=0
					count_3=0
					count_4=0
				if s%100000==0:
					i+=1
					pylab.figure(i)		#Plots total energy
					pylab.plot(E_, marker='.', linestyle='None')
					#pylab.ylim(555000,565000)
					pylab.show()
			
		return T,d

	A=time()
	T,d=run()
	'''for i in range(len(T[:])):
		for l in range(len(T[:][i])):
			x,y= T[:][i][l],  d_make(d,i)[l]
			if x!='N/A' and y!='N/A':
				print x,y'''
	print time()-A
