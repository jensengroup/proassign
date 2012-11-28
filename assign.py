import random
from numpy import exp, sqrt, pi, log, average, std, zeros, savetxt
import pylab
#import copy
import reader
from time import time, clock


steps=20000000 #maxsteps to be run

Names=['1CEX']	#['1A2P', '1CEX', '1DMB','1GZI', '1HCB', '1UBQ', '2GB1'] #Name used for plot and datalocation. Real data.
#Names=['1A2P','1CEX','1DMB','1GZI','1HCB','1LDE', '1LLW','1RLC','1UBQ','2GB1','2YP2','2ZKR','3S9I','3TTQ','3VUO','4AMC','4D94','4HOL'] #Synthetic data cs.
CS_pairs=['H', 'HA', 'N', 'C', 'CA', 'CB']	#set of CS pairs that must exist in both datasets. ('HA' refers to HA3 on glycine. HA2 can be added for shiftx2 only.)
Note='real_NA' #text will be added to generated filenames

theory_input='shiftx2' # 'shiftx2' or 'camshift'

CS_sigma_cam={'H':0.56,'HA':0.28,'N':3.01,'C':1.38, 'CA':1.3, 'CB':1.36}		#dictionary with value=variance for key=CS type. #Kohlhoff 2009

CS_sigma_cs={'N':1.12, 'CA':0.44, 'CB':0.52, 'C':0.98, 'HN':0.17, 'HA':0.12, 'H':0.17, 'CD':0.98, 'CG': 0.52, 'HB1':0.12, 'HB2':0.12, 'HB3':0.12} 	#shiftx2, CD and CG variance set to same as similar shift types, Han, B  Liu, Y, 2011


#######################################################################



if theory_input=='camshift':		# Determines which sigma dictionary to use
	CS_sigma=CS_sigma_cam
elif theory_input=='shiftx2':
	CS_sigma=CS_sigma_cs



def getE(t,d_e): #calculates total Energy difference of the two systems. 
	E=0
	for i in range(len(t)):
		E+=getdE_part(t,i,d_e)
	return E


def getdE_part(t,i,d_e):		#get the current energy difference of t and d_e at site i.
	dE=0
	for n in range(len(CS_pairs)):
		try:
			d=t[i][1][n]-d_e[i][1][n]
			dE += (d**2) /(2 * CS_sigma[CS_pairs[n]]**2)
		except TypeError:	#pass if there's one or both of the probed chemical shifts that's missing.
			pass
	return dE


def getdE(t,i,j,d_e):		#get energy difference for interchanging i and j
	dE=0
	for n in range(len(CS_pairs)):
		try:
			d_ij=t[i][1][n]-d_e[j][1][n]
			d_ji=t[j][1][n]-d_e[i][1][n]
			d_ii=t[i][1][n]-d_e[i][1][n]
			d_jj=t[j][1][n]-d_e[j][1][n]
			dE += ((d_ij)**2+(d_ji)**2-(d_ii)**2-(d_jj)**2) /(2 * CS_sigma[CS_pairs[n]]**2)
		except TypeError:
			pass
	return dE

def pair(t,d_e,d_t,d_le):		#t = a set of CS types to be used
	d__e=dict()	
	d__t=dict()
	for i in range(d_le):
		temp_e=[]
		temp_t=[]
		for l in t:
			try:
				temp_e.append(d_e[(i, l)])
			except KeyError:
				temp_e.append('N/A')
			temp_t.append(d_t[(i, l)])

		d__e[i]=temp_e
		d__t[i]=temp_t

	d__e=d__e.items()
	d__e.sort()
	d__t=d__t.items()
	d__t.sort()
	return d__e, d__t

def histogram(d_t,A,d_le, count): #The most probable assignment of the first residue will be the maximum value of row 0 (A[0,:]).
	for l in range(d_le):
		A[l,d_t[l][0]]+=count
	return A

def run():
	i=0
	for n in Names:
		energy=open(n+'_'+theory_input+'_'+Note+'_energy.txt', 'w')
		if theory_input=='camshift':
			d_t, d_le = reader.cam_(n+'.cam')	#camshift filename
			CS_sigma=CS_sigma_cam
		elif theory_input=='shiftx2':
			d_t,d_le = reader.cs_(n+'.cs') # name of shiftx2 file
			CS_sigma=CS_sigma_cs
		else:
			print 'Error in theory_input'
			break

		d_e, structure, d_le = reader.str_(n+'.str') # name of srt file
		#d_e, d_le=reader.cs_(n+'.cs') #use if using synthetic data
		d_e, d_t = pair(CS_pairs,d_e, d_t,d_le)
		'''for d in d_t:			#Randomizes data with the default variances (only for synthetic data)
			for x in range(len(CS_pairs)):
				try:
					d[1][x]=random.gauss(d[1][x],CS_sigma[CS_pairs[x]])
				except TypeError:
					pass'''
		A=zeros((d_le,d_le))
		flag=False
		E_=[]
		E_ini=getE(d_t,d_e)	
		random.shuffle(d_t)
		percent=0
		E_tot=getE(d_t,d_e)-E_ini
		s=0
		S=0
		count=0	#count for how many times the same configuration occurs before a change.
		while s<steps:
			s+=1
			S+=1
			i, j = random.randint(0,len(d_t)-1), random.randint(0,len(d_t)-1)
			dE=getdE(d_t,i,j,d_e)
			if dE<=0:
				if flag==True:
					A=histogram(d_t,A,d_le,count)
				E_before=getdE_part(d_t,i,d_e)+getdE_part(d_t,j,d_e)
				d_t[i], d_t[j] = d_t[j], d_t[i]
				E_after=getdE_part(d_t,i,d_e)+getdE_part(d_t,j,d_e)
				E_tot+=E_after-E_before
				count=1
			else:
				p=random.random()
				if p<exp(-dE):
					if flag==True:
						A=histogram(d_t,A,d_le,count)
					E_before=getdE_part(d_t,i,d_e)+getdE_part(d_t,j,d_e)
					d_t[i], d_t[j] = d_t[j], d_t[i]
					E_after=getdE_part(d_t,i,d_e)+getdE_part(d_t,j,d_e)
					E_tot+=E_after-E_before
					count=1
			E_.append(E_tot)
			count +=1
			if flag!=True:				
				if s== 8*d_le**2:		# number of steps that you're sure the energy stabilize prior to.
					s=max(steps-s,1)	# run s steps after stabilization (prior step)
					flag=True

			energy.write(str(E_tot)+'\n')
			if S%5000==0:
				print S, n
	savetxt(n+'_'+theory_input+'_'+Note+'.txt',A)
	energy.close()
	i+=1
	pylab.figure(i)		#Plots total energy
	pylab.plot(E_, marker='.', linestyle='None')
	pylab.xlabel('steps')
	pylab.ylabel('Energy difference')
	pylab.ylim(min(E_)-abs(min(E_))*1.5,abs(max(E_[100*d_le:]))*1.5)
	pylab.savefig(n+'_'+theory_input+Note+'_'+'_energy')

A=time()
run()
print time()-A
