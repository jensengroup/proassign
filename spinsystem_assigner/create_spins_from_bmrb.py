def HSQC(d,CS_pairs):
	t=[]
	for i in range(len(CS_pairs)):
		if CS_pairs[i]=='N':
			for l in range(len(CS_pairs)):
				if CS_pairs[l]=='H' or CS_pairs[l]=='HA2' or CS_pairs[l]=='HA3':
					for n in d:
						if n[i]!='N/A' and n[l]!='N/A':
							t.append([n[i],n[l]])
	data=open('HSQC.list','wb')
	data.write('\n')
	data.write('\n')
	for i in t:
		data.write('?-? %s %s \n'% (i[0],i[1]))

def HNCA(d,CS_pairs):
	t=[]
	for i in range(len(CS_pairs)):
		if CS_pairs[i]=='N':
			for l in range(len(CS_pairs)):
				if CS_pairs[l]=='H':
					for m in range(len(CS_pairs)):
						if CS_pairs[m]=='CA':
							for n in range(len(d)):
								try:
									t.append([float(d[n][m]),float(d[n][i]),float(d[n][l])])
								except ValueError:
									pass
								try:
									if n!=0:
										t.append([float(d[n-1][m]),float(d[n][i]),float(d[n][l])])
								except ValueError:
									pass
	data=open('HNCA.list','wb')
	data.write('\n')
	data.write('\n')
	for i in t:
		data.write('?-? %s %s %s\n'% (i[0],i[1], i[2]))

def HNCACB(d,CS_pairs):
	t=[]
	for i in range(len(CS_pairs)):
		if CS_pairs[i]=='N':
			for l in range(len(CS_pairs)):
				if CS_pairs[l]=='H':
					for m in range(len(CS_pairs)):
						if CS_pairs[m]=='CA' or CS_pairs[m]=='CB':
							for n in range(len(d)):
								try:
									t.append([float(d[n][m]),float(d[n][i]),float(d[n][l])])
								except ValueError:
									pass
								try:
									if n!=0:
										t.append([float(d[n-1][m]),float(d[n][i]),float(d[n][l])])
								except ValueError:
									pass
	data=open('HNCACB.list','wb')
	data.write('\n')
	data.write('\n')
	for i in t:
		data.write('?-? %s %s %s\n'% (i[0],i[1], i[2]))

def HNCO(d,CS_pairs):
	t=[]
	for i in range(len(CS_pairs)):
		if CS_pairs[i]=='N':
			for l in range(len(CS_pairs)):
				if CS_pairs[l]=='H':
					for m in range(len(CS_pairs)):
						if CS_pairs[m]=='C':
							for n in range(len(d)):
								try:
									if n!=0:
										t.append([float(d[n-1][m]),float(d[n][i]),float(d[n][l])])
								except ValueError:
									pass
	data=open('HNCO.list','wb')
	data.write('\n')
	data.write('\n')
	for i in t:
		data.write('?-? %s %s %s\n'% (i[0],i[1], i[2]))

def HNcoCA(d,CS_pairs):
	t=[]
	for i in range(len(CS_pairs)):
		if CS_pairs[i]=='N':
			for l in range(len(CS_pairs)):
				if CS_pairs[l]=='H':
					for m in range(len(CS_pairs)):
						if CS_pairs[m]=='CA':
							for n in range(len(d)):
								try:
									if n!=0:
										t.append([float(d[n-1][m]),float(d[n][i]),float(d[n][l])])
								except ValueError:
									pass
	data=open('HNcoCA.list','wb')
	data.write('\n')
	data.write('\n')
	for i in t:
		data.write('?-? %s %s %s\n'% (i[0],i[1], i[2]))

def HNcaCACB(d,CS_pairs):
	t=[]
	for i in range(len(CS_pairs)):
		if CS_pairs[i]=='N':
			for l in range(len(CS_pairs)):
				if CS_pairs[l]=='H':
					for m in range(len(CS_pairs)):
						if CS_pairs[m]=='CA' or CS_pairs[m]=='CB':
							for n in range(len(d)):
								try:
									t.append([float(d[n][m]),float(d[n][i]),float(d[n][l])])
								except ValueError:
									pass
								try:
									if n!=0:
										t.append([float(d[n-1][m]),float(d[n][i]),float(d[n][l])])
								except ValueError:
									pass
	data=open('HNcaCACB.list','wb')
	data.write('\n')
	data.write('\n')
	for i in t:
		data.write('?-? %s %s %s\n'% (i[0],i[1], i[2]))

def HNcoHA(d,CS_pairs):
	t=[]
	for i in range(len(CS_pairs)):
		if CS_pairs[i]=='N':
			for l in range(len(CS_pairs)):
				if CS_pairs[l]=='H':
					for m in range(len(CS_pairs)):
						if CS_pairs[m]=='HA' or CS_pairs[m]=='HA2' or CS_pairs[m]=='HA3':
							for n in range(len(d)):
								try:
									if n!=0:
										t.append([float(d[n-1][m]),float(d[n][i]),float(d[n][l])])
								except ValueError:
									pass
	data=open('HNcoHA.list','wb')
	data.write('\n')
	data.write('\n')
	for i in t:
		data.write('?-? %s %s %s\n'% (i[0],i[1], i[2]))

def HNHA(d,CS_pairs):
	t=[]
	for i in range(len(CS_pairs)):
		if CS_pairs[i]=='N':
			for l in range(len(CS_pairs)):
				if CS_pairs[l]=='H':
					for n in range(len(d)):
						for m in range(len(CS_pairs)):
							if CS_pairs[m]=='HA' or CS_pairs[m]=='HA2' or CS_pairs[m]=='HA3':
								try:
									t.append([float(d[n][m]),float(d[n][i]),float(d[n][l])])
								except ValueError:
									pass
								try:
									if n!=0:
										t.append([float(d[n-1][m]),float(d[n][i]),float(d[n][l])])
								except ValueError:
									pass
	data=open('HNHA.list','wb')
	data.write('\n')
	data.write('\n')
	for i in t:
		data.write('?-? %s %s %s\n'% (i[0],i[1], i[2]))

def run(d,CS_pairs):
	HSQC(d,CS_pairs)
	HNCA(d,CS_pairs)
	HNCACB(d,CS_pairs)
	HNCO(d,CS_pairs)
	HNcoCA(d,CS_pairs)
	HNcaCACB(d,CS_pairs)
	HNcoHA(d,CS_pairs)
	HNHA(d,CS_pairs)
