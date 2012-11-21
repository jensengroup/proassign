import random
from numpy import exp, sqrt, pi, log, average, std, zeros, savetxt
import pylab
#import copy
import reader_array as reader
from time import time, clock


steps=5000 #maxsteps to be run
repeat=1
#Names=['1CEX']#'1GZI', '1A2P', '1CEX', '1HCB', '1DMB', '1ubq', '2gb1'] #Name used for plot and datalocation.
Names=['1LLW']#'1LDE', '1LLW','1RLC','2YP2','2ZKR','3S9I','3TTQ','3VUO','4AMC','4D94','4HOL'] #Name used for plot and datalocation.
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


def getdE_part(t,i,d_e):		#get the current energy difference of t and d_e at site i.
	dE=0
	for n in range(len(CS_pairs)):
		if t[i][n+1]>0 and d_e[i][n+1]>0:
			d=t[i][n+1]-d_e[i][n+1]
			dE += (d**2) /(2 * CS_sigma[CS_pairs[n]]**2)
		else:
			pass
	return dE


def getdE(t,i,j,d_e):		#get energy difference for interchanging i and j
	dE=0
	for n in range(len(CS_pairs)):
		if t[i][n+1]>0 and t[j][n+1]>0 and d_e[i][n+1]>0 and d_e[j][n+1]>0:
			d_ij=t[i][n+1]-d_e[j][n+1]
			d_ji=t[j][n+1]-d_e[i][n+1]
			d_ii=t[i][n+1]-d_e[i][n+1]
			d_jj=t[j][n+1]-d_e[j][n+1]
			dE += ((d_ij)**2+(d_ji)**2-(d_ii)**2-(d_jj)**2) /(2 * CS_sigma[CS_pairs[n]]**2)
		else: 
			pass
	return dE

def pair(t,d_e,d_t,d_le):		#t = a set of CS types to be used
	d__t=zeros((d_le,1+len(CS_pairs)))
	d__e=zeros((d_le,1+len(CS_pairs)))
	for i in range(d_le):
		temp_e=[]
		temp_t=[]
		for l in range(len(t)):
			try:
				d__e[i][l+1]=float(d_e[(i, t[l])])
			except KeyError:
				d__e[i][l+1]=-1
			try:
				d__t[i][l+1]=float(d_t[(i, t[l])])
			except (KeyError, ValueError):
				try:
					l=-1
					d__t[i][l+1]=float(d_t[(i, t[l])])
				except (KeyError, ValueError):
					d__t[i][l+1]=1

	return d__e, d__t

def histogram(d_t,A,d_le): #slow. Can do faster but will use more hdd space and infinite ram to parse.
	for l in range(d_le):
		A[l,d_t[l][0]]+=1
	return A

def run():
	i=0
	for n in Names:
		#energy=open(n+'_'+theory_input+'_energy.txt', 'w')
		if theory_input=='camshift':
			d_t = reader.cam_(n+'.cam')	#camshift filename
			CS_sigma=CS_sigma_cam
		elif theory_input=='shiftx2':
			d_t, d_le = reader.cs_(n+'.cs') # name of shiftx2 file
			CS_sigma=CS_sigma_cs
		else:
			print 'Error in theory_input'
			break

		#d_e, structure, d_le = reader.str_(n+'.str') # name of srt file
		d_e, d_le=reader.cs_(n+'.cs')

		d_e, d_t = pair(CS_pairs,d_e, d_t,d_le)
		for d in d_t:			#Randomizes data with the default variances
			for x in range(len(CS_pairs)):
				try:
					d[x+1]=random.gauss(d[x+1],CS_sigma[CS_pairs[x]])
				except TypeError:
					pass
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
				if s==1:
					flag=True
				s+=1
				S+=1
				i, j = random.randint(0,len(d_t)-1), random.randint(0,len(d_t)-1)
				dE=getdE(d_t,i,j,d_e)
				if dE<=0:
					E_before=getdE_part(d_t,i,d_e)+getdE_part(d_t,j,d_e)
					d_t[i][1:], d_t[j][1:] = d_t[j][1:], d_t[i][1:]
					E_after=getdE_part(d_t,i,d_e)+getdE_part(d_t,j,d_e)
					E_tot+=E_after-E_before
				else:
					p=random.random()
					if p<exp(-dE):
						E_before=getdE_part(d_t,i,d_e)+getdE_part(d_t,j,d_e)
						d_t[i][1:], d_t[j][1:] = d_t[j][1:], d_t[i][1:]
						E_after=getdE_part(d_t,i,d_e)+getdE_part(d_t,j,d_e)
						E_tot+=E_after-E_before
				'''if 20*s%steps == 0:
					percent+=5
					print percent,"%"'''
				E_.append(E_tot)
				if flag!=True:				
					if s%10000==0:
						if float(E_[s-1]/E_[1*s/2])>0.99:
							s=max(steps-s,1)
							flag=True
				if flag==True:
					A=histogram(d_t,A,d_le)
				#energy.write(str(E_tot)+'\n')
				if S%5000==0:
					print S, n
		savetxt(n+'__'+theory_input+'.txt',A)
		#energy.close()
		i+=1
		pylab.figure(i)		#Plots total energy
		pylab.plot(E_, marker='.', linestyle='None')
		pylab.xlabel('steps')
		pylab.ylabel('Energy difference')
		pylab.ylim(min(E_)-abs(E_[-1]),abs(E_[-1]*10))
		pylab.savefig('Assignment %s' %n)

A=time()
run()
print time()-A #10.3 5k
