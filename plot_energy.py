import pylab

#energy plots for the energy dump from assign.py

def lines(filename):				#Returns list_[lines][elements] from the source data.
	data = open(filename, 'r')
	list_ = []
	for line in data:
		list_.append(float(line.split()[0]))
	data.close()
	return list_


A=['1CEX']#[ '1A2P', '1CEX', '1DMB','1GZI', '1HCB', '1UBQ', '2GB1']

for i in range(len(A)):
	E_= lines(A[i]+'_shiftx2_real_NA_energy.txt')
	pylab.figure(i)		#Plots total energy
	pylab.plot(E_[:100000], marker='.', linestyle='None')
	pylab.xlabel('steps')
	pylab.ylabel('Energy difference')
	#pylab.ylim(-1000,2000)
	pylab.show(A[i]+'_energy')
