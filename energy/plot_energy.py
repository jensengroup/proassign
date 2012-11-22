import pylab

def lines(filename):				#Returns list_[lines][elements] from the source data.
	data = open(filename, 'r')
	list_ = []
	for line in data:
		list_.append(float(line.split()[0]))
	data.close()
	return list_

E_= lines('2gb1_shiftx2_shuffle_energy.txt')

pylab.figure(1)		#Plots total energy
pylab.plot(E_[:10000], marker='.', linestyle='None')
pylab.xlabel('steps')
pylab.ylabel('Energy difference')
pylab.ylim(-1000,2000)
pylab.savefig('energy')
