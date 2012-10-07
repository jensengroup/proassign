import math
import numpy
import string
import sys
import os

ALPHABET_UPPER_CASE = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
ALPHABET_LOWER_CASE = "abcdefghijklmnoqrstuvwxyz"
ALPHABET_NUMBERS = "0123456789"

def one_to_three(letter):
	if letter == "A": return "ALA"
	elif letter == "G": return "GLY"
	elif letter == "T": return "THR"
	elif letter == "S": return "SER"
	elif letter == "C": return "CYS"
	elif letter == "V": return "VAL"
	elif letter == "L": return "LEU"
	elif letter == "I": return "ILE"
	elif letter == "M": return "MET"
	elif letter == "P": return "PRO"
	elif letter == "F": return "PHE"
	elif letter == "Y": return "TYR"
	elif letter == "D": return "ASP"
	elif letter == "E": return "GLU"
	elif letter == "W": return "TRP"
	elif letter == "N": return "ASN"
	elif letter == "Q": return "GLN"
	elif letter == "H": return "HIS"
	elif letter == "R": return "ARG"
	else: return "UNK"


def rmsd(data1, data2, atom_type):
	rmsd = 0.0
	n_count = 0.0
	for atom1 in data1:
		if atom1[2] == atom_type:
			for atom2 in data2:
				if (atom1[0] == atom2[0]) and (atom1[1] == atom2[1]) and (atom1[2] == atom2[2]):
					
#					print atom2
					n_count = n_count + 1.0
					rmsd = rmsd + (atom1[3] - atom2[3])*(atom1[3] - atom2[3])
#					print atom1, atom2, atom1[3] - atom2[3], n_count, rmsd/n_count
	rmsd = math.sqrt(rmsd/n_count)
#	print n_count, 
	if n_count == 0:
		return 999.9
	else:
		return rmsd

def carbon_this_reftype_to_dss(molecule, reftype):
	if reftype == "DSS": offset = 0
	elif reftype == "TMS": offset = - 1.7
	elif reftype == "TSP_PH7": offset = -0.12
	elif reftype == "TSP_PH3": offset = -0.22
	elif reftype == "NaAc": offset = 26.1
	elif reftype == "Acetone": offset = 33.0
	elif reftype == "Dioxane": offset = 69.3
	else: 
		print "Unknown carbon reference: ", reftype, ". Aborting ..."
		sys.exit()
	for atom in molecule:
		if "C" in atom[2]: atom[3] = atom[3] - offset
	return

def hydrogen_this_reftype_to_dss(molecule, reftype):
	if reftype ==   "DSS":     offset = 0
	if reftype ==   "TMS":     offset = 0
	elif reftype == "TSP_PH7": offset = -0.015
	elif reftype == "TSP_PH3": offset = 0.003
	elif reftype == "Acetone": offset = 2.218
	elif reftype == "Dioxane": offset = 3.750
	elif reftype == "HDO":     offset = 4.766
	else: 
		print "Unknown hydrogen reference: ", reftype, ". Aborting ..."
		sys.exit()
	for atom in molecule:
		if "H" in atom[2]: atom[3] = atom[3] - offset
	return

def nitrogen_this_parallel_reftype_to_nh3(molecule, reftype):
	if reftype ==   "NH3":         offset = 0
	elif reftype == "NH4NO3":      offset = 21.0
	elif reftype == "NH4Cl":       offset = 23.6
	elif reftype == "Urea":        offset = 77.0
	elif reftype == "CH3NO2_1to1": offset = 379.8
	elif reftype == "CH3NO2_neat": offset = 381.7
	else: 
		print "Unknown nitrogen reference: ", reftype, ". Aborting ..."
		sys.exit()
	for atom in molecule:
		if "N" in atom[2]: atom[3] = atom[3] - offset
	return

def nitrogen_this_perpendicular_reftype_to_nh3(molecule, reftype):
	if reftype ==   "NH3":         offset = 0
	elif reftype == "NH4NO3":      offset = 21.6
	elif reftype == "NH4Cl":       offset = 24.9
	elif reftype == "Urea":        offset = 76.3
	elif reftype == "CH3NO2_1to1": offset = 379.6
	elif reftype == "CH3NO2_neat": offset = 380.2
	else: 
		print "Unknown nitrogen reference: ", reftype, ". Aborting ..."
		sys.exit()
	for atom in molecule:
		if "N" in atom[2]: atom[3] = atom[3] - offset
	return

def nitrogen_this_magicangle_reftype_to_nh3(molecule, reftype):
	if reftype ==   "NH3":         offset = 0
	elif reftype == "NH4NO3":      offset = 22.8
	elif reftype == "NH4Cl":
		print "Sorry, no data for a magic angle reference of the type: ", reftype, ". Aborting ..."
		sys.exit()
	elif reftype == "Urea":        offset = 79.4
	elif reftype == "CH3NO2_1to1": offset = 378.1
	elif reftype == "CH3NO2_neat": offset = 381.9
	else: 
		print "Unknown nitrogen reference: ", reftype, ". Aborting ..."
		sys.exit()
	for atom in molecule:
		if "N" in atom[2]: atom[3] = atom[3] - offset
	return


def readfile(file_name, file_type):
	if not os.path.exists(file_name):
		print "ERROR: Could not open file", file_name
		sys.exit()
#	else:
#		print "Found file", file_name
	if file_type == "brmb":
#		rpint "File type is BRMB"
		return readfile_brmb(file_name)
	elif file_type == "shiftx":	
#		print "File type is SHIFTX"
		return readfile_shiftx(file_name)
	elif file_type == "shiftx2":	
#		print "File type is .cs SHIFTX2"
		return readfile_shiftx2(file_name)
	elif file_type == "shifty":	
#		print "File type is .shifty SHIFTY (SHIFTX2)"
		return readfile_shifty(file_name)
	elif file_type == "sparta":	
#		print "File type is SPARTA std. output"
		return readfile_sparta(file_name)
	elif file_type == "spartaplus":	
#		print "File type is SPARTA+ std. output"
		return readfile_spartaplus(file_name)
	elif file_type == "camshift":	
#		print "File type is CamShift std. output"
		return readfile_camshift(file_name)
	elif file_type == "qdb":	
#		print "File type is .qdb SHIFTS"
		return readfile_qdb(file_name)
	elif file_type == "emp":	
#		print "File type is .emp SHIFTS"
		return readfile_emp(file_name)
	elif file_type == "cheshift":	
#		print "File type is .out CheSHIFT"
		return readfile_cheshift(file_name)
	elif file_type == "proshift":	
#		print "File type is .pro PROSHIFT"
		return readfile_proshift(file_name)
	elif file_type == "padawan":	
#		print "File type is verbose Padawan std. output"
		return readfile_padawan(file_name)
	else:
		print "ERROR: Unknown filetype", file_type


def readfile_padawan(file_name):
	atoms_list = []
        file_text = open(file_name)
        file_line = file_text.readlines()
	mode = "start"
	read_counter = 0
	for line in file_line:
		if mode == "read_cs":
			if (len(line) > 5) and ("=" not in line):
#				print line
				number = int(line[3]+line[4]+line[5])
				name = line[0]+line[1]+line[2]
				cs = float(string.split(line)[2])
				atom_data = [number, name, "H", cs]
				atoms_list.append(atom_data)
#				print atom_data, "\n"
		if mode == "start":
			if "Residue:    Experimental:[ppm]      Calculated:[ppm]" in line:
				mode = "read_cs"
#	print atoms_list
	return atoms_list


def readfile_proshift(file_name):
	atoms_list = []
        file_text = open(file_name)
        file_line = file_text.readlines()
	mode = "start"
	read_counter = 0
	for line in file_line:
		if (len(line) > 5):
			if string.split(line)[0] == "SHIFT":
#				print line,
				atom_data = [int(string.split(line)[5]), string.split(line)[3], string.split(line)[2], float(string.split(line)[6])]
#				print atom_data, "\n"
				atoms_list.append(atom_data)
#	print atoms_list
	return atoms_list



def readfile_cheshift(file_name):
	atoms_list = []
        file_text = open(file_name)
        file_line = file_text.readlines()
	mode = "start"
	read_counter = 0
	for line in file_line:
		if mode == "read_cs":
			if (len(line) > 5) and ("Date" not in line):
#				print line
				atom_data = [int(string.split(line)[0]), string.split(line)[1], "CA", float(string.split(line)[2])]
				if (atom_data[3] > 0) and (atom_data[3]< 998.0): atoms_list.append(atom_data)
			else:
				mode = "stop"
#			print atom_data, "\n"
		if mode == "start":
			if "Chemical shifts" in line:
				mode = "read_cs"
#	print atoms_list
	carbon_this_reftype_to_dss(atoms_list, "TMS")
	return atoms_list

def readfile_emp(file_name):
	atoms_list = []
        file_text = open(file_name)
        file_line = file_text.readlines()
	mode = "start"
	read_counter = 0
	for line in file_line:
		if mode == "read_cs":
			if (len(line) > 5):
#				print line
				atom_data = [int(string.split(line)[2]), string.split(line)[1], string.split(line)[3], float(string.split(line)[8])]
				atoms_list.append(atom_data)
#			print atom_data, "\n"
		if mode == "start":
			if "=============================================================================" in line:
				read_counter = read_counter + 1
			if read_counter == 4:
				mode = "read_cs"
#	print atoms_list
	return atoms_list


def readfile_qdb(file_name):
	atoms_list = []
        file_text = open(file_name)
        file_line = file_text.readlines()
	for line in file_line:
		if (len(line) > 5) and ( "#" not in line):
#q			print line,
			atom_data = [int(string.split(line)[3]), string.split(line)[2], string.split(line)[4], float(string.split(line)[13])]
			if atom_data[2] == "Ca": atom_data[2] = "CA"
			if atom_data[2] == "Cb": atom_data[2] = "CB"
			atoms_list.append(atom_data)
#			print atom_data, "\n"
	return atoms_list


def readfile_camshift(file_name):
	atoms_list = []
        file_text = open(file_name)
        file_line = file_text.readlines()
	mode = "wait"
	for line in file_line:
		if mode == "read":
			if len(line) > 5:
#				print line,
				atom_data = [int(string.split(line)[0]), string.split(line)[1],
					     "HA", float(string.split(line)[2])]
				if atom_data[3] > 0.0: atoms_list.append(atom_data)

				atom_data = [int(string.split(line)[0]), string.split(line)[1],
					     "CA", float(string.split(line)[3])]
				if atom_data[3] > 0.0: atoms_list.append(atom_data)

				atom_data = [int(string.split(line)[0]), string.split(line)[1],
					     "H", float(string.split(line)[4])]
				if atom_data[3] > 0.0: atoms_list.append(atom_data)

				atom_data = [int(string.split(line)[0]), string.split(line)[1],
					     "N", float(string.split(line)[5])]
				if atom_data[3] > 0.0: atoms_list.append(atom_data)

				atom_data = [int(string.split(line)[0]), string.split(line)[1],
					     "C", float(string.split(line)[6])]
				if atom_data[3] > 0.0: atoms_list.append(atom_data)

				atom_data = [int(string.split(line)[0]), string.split(line)[1],
					     "CB", float(string.split(line)[7])]
				if atom_data[3] > 0.0: atoms_list.append(atom_data)
			else: mode = "stop"
		if mode == "wait":
			if "---------------------------------------------------------------------" in line:
				mode = "read"
#	print atoms_list
	return atoms_list

def readfile_spartaplus(file_name):
	atoms_list = []
        file_text = open(file_name)
        file_line = file_text.readlines()
	mode = "wait"
	for line in file_line:
		if mode == "read":
			if len(line) > 5:
				atom_data = [int(string.split(line)[0]), one_to_three(string.split(line)[1]),
					     string.split(line)[2], float(string.split(line)[4])]
#				print atom_data
#				print line
				if atom_data[2] == "HN": atom_data[2] = "H"
				atoms_list.append(atom_data)
		if mode == "wait":
			if "FORMAT %4d %4s %4s %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f" in line:
				mode = "read"
#	print atoms_list
	return atoms_list

def readfile_spartaplus(file_name):
	atoms_list = []
        file_text = open(file_name)
        file_line = file_text.readlines()
	mode = "wait"
	for line in file_line:
		if mode == "read":
			if len(line) > 5:
				atom_data = [int(string.split(line)[0]), one_to_three(string.split(line)[1]),
					     string.split(line)[2], float(string.split(line)[5])]
#				print atom_data
#				print line
				if atom_data[2] == "HN": atom_data[2] = "H"
				atoms_list.append(atom_data)
		if mode == "wait":
			if "FORMAT %4d %4s %4s %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f" in line:
				mode = "read"
#	print atoms_list
	return atoms_list

def readfile_sparta(file_name):
	atoms_list = []
        file_text = open(file_name)
        file_line = file_text.readlines()
	mode = "wait"
	for line in file_line:
		if mode == "read":
			if len(line) > 5:
#				print line,
				atom_data = [int(string.split(line)[6]), one_to_three(string.split(line)[7]),
					     "N", float(string.split(line)[0])]
				if atom_data[3] < 9998.00: atoms_list.append(atom_data)

				atom_data = [int(string.split(line)[6]), one_to_three(string.split(line)[7]),
					     "HA", float(string.split(line)[1])]
				if atom_data[3] < 9998.00: atoms_list.append(atom_data)

				atom_data = [int(string.split(line)[6]), one_to_three(string.split(line)[7]),
					     "C", float(string.split(line)[2])]
				if atom_data[3] < 9998.00: atoms_list.append(atom_data)

				atom_data = [int(string.split(line)[6]), one_to_three(string.split(line)[7]),
					     "CA", float(string.split(line)[3])]
				if atom_data[3] < 9998.00: atoms_list.append(atom_data)

				atom_data = [int(string.split(line)[6]), one_to_three(string.split(line)[7]),
					     "CB", float(string.split(line)[4])]
				if atom_data[3] < 9998.00: atoms_list.append(atom_data)
				atom_data = [int(string.split(line)[6]), one_to_three(string.split(line)[7]),
					     "H", float(string.split(line)[5])]
				if atom_data[3] < 9998.00: atoms_list.append(atom_data)
			else: mode = "stop"
		if mode == "wait":
			if "N       HA        C       CA       CB        H" in line:
				mode = "read"
#	print atoms_list
	return atoms_list


def readfile_shiftx2(file_name):
	atoms_list = []
        file_text = open(file_name)
        file_line = file_text.readlines()
	for line in file_line:
		if len(line) > 5:
#			print line,
			atom_data = [int(string.split(line,",")[0]), one_to_three(string.split(line,",")[1]), string.split(line,",")[2], float(string.split(line,",")[3])]
			atoms_list.append(atom_data)
#			print atom_data, "\n"
	return atoms_list


def readfile_shifty(file_name):
	atoms_list = []
        file_text = open(file_name)
        file_line = file_text.readlines()
	for line in file_line:
		if len(line) > 5:
#			print line,
			atom_data = [int(string.split(line,",")[0]), one_to_three(string.split(line,",")[1]), string.split(line,",")[2], float(string.split(line,",")[3])]
			atoms_list.append(atom_data)
#			print atom_data, "\n"
#	print atoms_list
	return atoms_list



def readfile_shiftx(file_name):
	atoms_list = []
        file_text = open(file_name)
        file_line = file_text.readlines()
	mode = "start"
	for line in file_line:
		if mode == "side_chain":
					atom_data = [ int(string.split(line)[0]), one_to_three(string.split(line)[1]), "HB", float(string.split(line)[4])]
					if string.split(line)[4] != "0.00": atoms_list.append(atom_data)
					atom_data = [ int(string.split(line)[0]), one_to_three(string.split(line)[1]), "HB2", float(string.split(line)[5])]
					if string.split(line)[5] != "0.00": atoms_list.append(atom_data)
					atom_data = [ int(string.split(line)[0]), one_to_three(string.split(line)[1]), "HB3", float(string.split(line)[6])]
					if string.split(line)[6] != "0.00": atoms_list.append(atom_data)
					atom_data = [ int(string.split(line)[0]), one_to_three(string.split(line)[1]), "HD1", float(string.split(line)[7])]
					if string.split(line)[7] != "0.00": atoms_list.append(atom_data)
					atom_data = [ int(string.split(line)[0]), one_to_three(string.split(line)[1]), "HD2", float(string.split(line)[8])]
					if string.split(line)[8] != "0.00": atoms_list.append(atom_data)
					atom_data = [ int(string.split(line)[0]), one_to_three(string.split(line)[1]), "HD21", float(string.split(line)[9])]
					if string.split(line)[9] != "0.00": atoms_list.append(atom_data)
					atom_data = [ int(string.split(line)[0]), one_to_three(string.split(line)[1]), "HD22", float(string.split(line)[10])]
					if string.split(line)[10] != "0.00": atoms_list.append(atom_data)
					atom_data = [ int(string.split(line)[0]), one_to_three(string.split(line)[1]), "HD3", float(string.split(line)[11])]
					if string.split(line)[11] != "0.00": atoms_list.append(atom_data)
					atom_data = [ int(string.split(line)[0]), one_to_three(string.split(line)[1]), "HE", float(string.split(line)[12])]
					if string.split(line)[12] != "0.00": atoms_list.append(atom_data)
					atom_data = [ int(string.split(line)[0]), one_to_three(string.split(line)[1]), "HE1", float(string.split(line)[13])]
					if string.split(line)[13] != "0.00": atoms_list.append(atom_data)
					atom_data = [ int(string.split(line)[0]), one_to_three(string.split(line)[1]), "HE2", float(string.split(line)[14])]
					if string.split(line)[14] != "0.00": atoms_list.append(atom_data)
					atom_data = [ int(string.split(line)[0]), one_to_three(string.split(line)[1]), "HE21", float(string.split(line)[15])]
					if string.split(line)[15] != "0.00": atoms_list.append(atom_data)
					atom_data = [ int(string.split(line)[0]), one_to_three(string.split(line)[1]), "HE22", float(string.split(line)[16])]
					if string.split(line)[16] != "0.00": atoms_list.append(atom_data)
					atom_data = [ int(string.split(line)[0]), one_to_three(string.split(line)[1]), "HE3", float(string.split(line)[17])]
					if string.split(line)[17] != "0.00": atoms_list.append(atom_data)
					atom_data = [ int(string.split(line)[0]), one_to_three(string.split(line)[1]), "HG", float(string.split(line)[18])]
					if string.split(line)[18] != "0.00": atoms_list.append(atom_data)
					atom_data = [ int(string.split(line)[0]), one_to_three(string.split(line)[1]), "HG1", float(string.split(line)[19])]
					if string.split(line)[19] != "0.00": atoms_list.append(atom_data)
					atom_data = [ int(string.split(line)[0]), one_to_three(string.split(line)[1]), "HG12", float(string.split(line)[20])]
					if string.split(line)[20] != "0.00": atoms_list.append(atom_data)
					atom_data = [ int(string.split(line)[0]), one_to_three(string.split(line)[1]), "HG13", float(string.split(line)[21])]
					if string.split(line)[21] != "0.00": atoms_list.append(atom_data)
					atom_data = [ int(string.split(line)[0]), one_to_three(string.split(line)[1]), "HG2", float(string.split(line)[22])]
					if string.split(line)[22] != "0.00": atoms_list.append(atom_data)
					atom_data = [ int(string.split(line)[0]), one_to_three(string.split(line)[1]), "HG3", float(string.split(line)[23])]
					if string.split(line)[23] != "0.00": atoms_list.append(atom_data)
					atom_data = [ int(string.split(line)[0]), one_to_three(string.split(line)[1]), "HZ", float(string.split(line)[24])]
					if string.split(line)[24] != "0.00": atoms_list.append(atom_data)
#					atom_data = [ int(string.split(line)[0]), one_to_three(string.split(line)[1]), #			print line,
		if mode == "backbone":
			if len(line) > 4:
				if "NUM RES  H    HA   HB   HB2  HB3  HD1  HD2  HD21 HD22 HD3  HE   HE1  HE2  HE21 HE22 HE3  HG   HG1  HG12 HG13 HG2  HG3  HZ" in line: mode = "side_chain"
				else:
#					print line,
					atom_data = [ int(string.split(line)[0]), one_to_three(string.split(line)[1]), "HA", float(string.split(line)[2])]
					if atom_data[3] > 1.0: atoms_list.append(atom_data)
					atom_data = [ int(string.split(line)[0]), one_to_three(string.split(line)[1]), "H", float(string.split(line)[3])]
					if atom_data[3] > 1.0: atoms_list.append(atom_data)
					atom_data = [ int(string.split(line)[0]), one_to_three(string.split(line)[1]), "N", float(string.split(line)[4])]
					if atom_data[3] > 1.0: atoms_list.append(atom_data)
					atom_data = [ int(string.split(line)[0]), one_to_three(string.split(line)[1]), "CA", float(string.split(line)[5])]
					if atom_data[3] > 1.0: atoms_list.append(atom_data)
					atom_data = [ int(string.split(line)[0]), one_to_three(string.split(line)[1]), "CB", float(string.split(line)[6])]
					if atom_data[3] > 1.0: atoms_list.append(atom_data)
					atom_data = [ int(string.split(line)[0]), one_to_three(string.split(line)[1]), "C", float(string.split(line)[7])]
					if atom_data[3] > 1.0: atoms_list.append(atom_data)
#					print atom_data
			
		if mode == "start":
			if "--- --- ------ ------ -------- ------- ------- --------" in line:
				mode = "backbone"
#	print atoms_list
	return atoms_list

def readfile_brmb(file_name):
	data_entries = ["_Atom_shift_assign_ID\n",
		        "_Residue_author_seq_code\n",
		        "_Residue_seq_code\n",
	        	"_Residue_label\n",
		        "_Atom_name\n",
		        "_Atom_type\n",
		        "_Chem_shift_value\n",
		        "_Chem_shift_value_error\n",
		        "_Chem_shift_ambiguity_code\n"]
	atoms_list = []
        file_text = open(file_name)
        file_line = file_text.readlines()
	data_structure = []
	lenght_data = 0
	n_count = 0
	mode = "loop"
	i_res_number = 0
	i_res_label = 0
	i_atom_label = 0
	i_chemical_shift = 0
	for line in file_line:
#SEARCH FOR loop_ STATEMENT
		if mode == "loop":
			if "loop_" in line:
				mode = "get"
#ASSIGN data_structure CONTENTS
		if mode == "get":
			if "_Atom_shift_assign_ID" in line:
				mode = "assign"
				n_count = 0
		if mode == "assign":
			if "_" in line:
				for i in data_entries:
					if "_Residue_seq_code\n" in line:
						i_res_number = n_count - 1 
					if "_Residue_label\n" in line:
						i_res_label = n_count - 1
					if "_Atom_name\n" in line:
						i_atom_label = n_count - 1
					if "_Chem_shift_value\n" in line:
						i_chemical_shift = n_count - 1
					if i in line:
						data_structure.append(i)
						n_count = n_count + 1
#				print line
			if len(line)< 10:
				mode = "read_cs"
				
		if mode == "read_cs":
#			print data_structure
#			print len(data_structure)
#			print i_res_number
#			print i_res_label
#			print i_atom_label
#			print i_chemical_shift
			if len(line) > 3:
				if "stop_" in line:
					break
				atom_data = [ int(string.split(line)[i_res_number]),string.split(line)[i_res_label],string.split(line)[i_atom_label],float(string.split(line)[i_chemical_shift])]
				atoms_list.append(atom_data)
#				print line,
	return atoms_list


