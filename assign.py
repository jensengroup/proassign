import random
from numpy import exp, sqrt, pi, log, average, std, zeros, savetxt
import pylab
#import copy
import reader
#from time import time, clock


steps=50000000
repeat=1
Names=['1GZI', '1A2P', '1CEX', '1HCB', '1DMB', '1ubq', '2gb1'] #Name used for plot and datalocation.
#Names=['1LDE', '1LLW','1RLC','2YP2','2ZKR','3S9I','3TTQ','3VUO','4AMC','4D94','4HOL'] #Name used for plot and datalocation.
CS_pairs=['H','HA','N',"C", 'CA', 'CB']			#set of CS pairs that must exist in both datasets. 


theory_input='shiftx2' # 'shiftx2' or 'camshift'

CS_sigma_cam={'H':0.6,'HA':0.35,'N':3.,'C':1.45, 'CA':1.35, 'CB':1.4}		#dictionary with value=variance for key=CS type. 
																			#{'H':0.6,'HA':0.35,'N':3.,"C'":1.45, 'CA':1.35, 'CB':1.4} 
																			#for camshift 30+ residue proteins Kohlhoff 2009
CS_sigma_cs={'N':1.12, 'CA':0.44, 'CB':0.52, 'C':0.98, 'HN':0.17, 'HA':0.12, 'H':0.17, 'CD':0.98, 'CG': 0.52, 'HB1':0.12, 'HB2':0.12, 'HB3':0.12} 	#shiftx2, CD and CG variance set to same as similar shift types, Han, B  Liu, Y, 2011
#######################################################################



if theory_input=='camshift':
	CS_sigma=CS_sigma_cam
elif theory_input=='shiftx2':
	CS_sigma=CS_sigma_cs



def getE(t,d_e): #calculates total Energy difference of the two systems. 
	E=0
	for i in range(len(t)):
		E+=getdE_part(t,i,d_e)
	return E



def is_small(i,j,n,value):
	if value[i][1][n] < 0.5 and value[j][1][n] < 0.5:
		return True
	return False 

def is_small_small(i,n,value):
	if value[i][1][n] < 0.5: return True
	return False
	
def getdE_part(t,i,d_e):		#get energy difference for interchanging i and j
	dE=0
	for n in range(len(CS_pairs)):
		#if t[i][1][n] !=0 and d_e[i][1][n] !=0:
		if not is_small_small(i,n,t) and not is_small_small(i,n,d_e):
			d=t[i][1][n]-d_e[i][1][n]
			dE += (d**2) /(2 * CS_sigma[CS_pairs[n]]**2)
	return dE


def getdE(t,i,j,d_e):		#get energy difference for interchanging i and j
	dE=0
	for n in range(len(CS_pairs)):
		#if t[i][1][n] > 1.0 and t[j][1][n] > 1.0 and d_e[i][1][n] > 1.0 and d_e[j][1][n] > 1.0:
		if not is_small(i,j,n,t) and not is_small(i,j,n,d_e): 
			d_ij=t[i][1][n]-d_e[j][1][n]
			d_ji=t[j][1][n]-d_e[i][1][n]
			d_ii=t[i][1][n]-d_e[i][1][n]
			d_jj=t[j][1][n]-d_e[j][1][n]
			dE += ((d_ij)**2+(d_ji)**2-(d_ii)**2-(d_jj)**2) /(2 * CS_sigma[CS_pairs[n]]**2)
	return dE

def pair(t,d_e,d_t,d_le):		#t = a set of CS types to be used
	d__e=dict()	
	d__t=dict()
	for i in range(d_le):
		temp_e=[]
		temp_t=[]
		for l in t:
			try:
				temp_e.append(float(d_e[(i, l)]))
			except KeyError:
				temp_e.append(0)
			try:
				temp_t.append(float(d_t[(i, l)]))
			except KeyError:
				try:
					l='HN'
					temp_t.append(float(d_t[(i, l)]))
				except KeyError:
					temp_t.append(0)

		d__e[i]=temp_e
		d__t[i]=temp_t

	d__e=d__e.items()
	d__e.sort()
	d__t=d__t.items()
	d__t.sort()
	return d__e, d__t

def histogram(d_t,A,d_le):
	for l in range(d_le):
		A[l,d_t[l][0]]+=1
	return A

def run():
	i=0
	for n in Names:
		#data=open(n+'_'+theory_input+'.txt', 'w')
		energy=open(n+'_'+theory_input+'_energy.txt', 'w')
		if theory_input=='camshift':
			d_t = reader.cam_(n+'.cam')	#camshift filename
			CS_sigma=CS_sigma_cam
		elif theory_input=='shiftx2':
			d_t = reader.cs_(n+'.cs') # name of shiftx2 file
			CS_sigma=CS_sigma_cs
		else:
			print 'Error in theory_input'
			break

		d_e, structure, d_le = reader.str_(n+'.str') # name of srt file
		#d_e=reader.cs_(n+'.cs')
		#d_e, d_t = pair(CS_pairs,d_e, d_t,len(d_t))
		
		'''for d in d_t:			#Randomizes data with the default variances
			for x in range(len(CS_pairs)):
				d[1][x]=random.gauss(d[1][x],CS_sigma[CS_pairs[x]])'''
		d_e, d_t = pair(CS_pairs,d_e, d_t,d_le)
		A=zeros((d_le,d_le))
		for r in range(repeat):
			flag=False
			E_=[]
			E_ini=getE(d_t,d_e)		
			random.shuffle(d_t)
			percent=0
			E_tot=getE(d_t,d_e)-E_ini+1e-10		#infinitisimal to prevent later division by zero.
			s=0
			S=0
			while s<steps:
				s+=1
				S+=1
				i, j = random.randint(0,len(d_t)-1), random.randint(0,len(d_t)-1)
				dE=getdE(d_t,i,j,d_e)
				if dE<=0:
					E_before=getdE_part(d_t,i,d_e)+getdE_part(d_t,j,d_e)
					d_t[i], d_t[j] = d_t[j], d_t[i]
					E_after=getdE_part(d_t,i,d_e)+getdE_part(d_t,j,d_e)
					E_tot+=E_after-E_before
				else:
					p=random.random()
					if p<exp(-dE):
						E_before=getdE_part(d_t,i,d_e)+getdE_part(d_t,j,d_e)
						d_t[i], d_t[j] = d_t[j], d_t[i]
						E_after=getdE_part(d_t,i,d_e)+getdE_part(d_t,j,d_e)
						E_tot+=E_after-E_before
				'''if 20*s%steps == 0:
					percent+=5
					print percent,"%"'''
				E_.append(E_tot)
				if flag!=True:				
					if s%1000==0:
						if float(E_[s-1]/E_[1*s/2])>0.99:
							s=max(steps-1*s,1)
							flag=True
				if flag==True:
					A=histogram(d_t,A,d_le)
				energy.write(str(E_tot)+'\n')
				if S%5000==0:
					print S, n
		savetxt(n+'_'+theory_input+'.txt',A)
		energy.close()
		i+=1
		
		pylab.figure(i)		#Plots total energy
		pylab.plot(E_, marker='.', linestyle='None')
		pylab.xlabel('steps')
		pylab.ylabel('Energy difference')
		pylab.ylim(min(E_)-abs(E_[-1]),abs(E_[-1]*10))
		pylab.savefig('Assignment %s' %n)


run()

