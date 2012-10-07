import numpy
from numpy.random import normal
from numpy.random import randint
from numpy.random import binomial

import src.multics as multics


def sample_nuisance(s, sp, ss, ssp, dlog):
    i = randint(0, 2)

    return s, sp, ss, ssp, dlog


def calc_energy():
    return 0


def calc_log_bias():
    return 0

def assemble_HN_matrix(data, length):
    b = numpy.zeros((length, 2))
    for atom in camshift_data:
        index = atom[0] - 1
        if atom[2] == 'H':
            b[index][0] = atom[3]
        if atom[2] == 'N':
            b[index][1] = atom[3]
    return b



## START OF PROGRAM

s    = [0.4, 2.9]
sp   = [0.4, 2.9]

ss   = [1.0, 1.0]
ssp  = [1.0, 1.0]

dlog = [1.0, 1.0]


n_steps = int(1E6)

# ref_data        = multics.readfile("bmr4317.str.corr", "brmb")
camshift_data   = multics.readfile("rn1ail.out.1.pdb.camshift", "camshift")

experimental = assemble_HN_matrix(camshift_data, 70)


## START OF ITERATIONS

for i in range(n_steps):

    sample_nuisance = binomial(n=1, p=.3)

    if sample_nuisance:
        s, sp, ss, ssp, dlog = sample_nuisance(s, sp, ss, ssp, dlog)

    else:
        a = randint(0, len(matrix))
        b = randint(0, len(matrix))
        pred[a], pred[b] = pred[a], pred[b]


    energy  = calc_energy()
    bias    = calc_log_bias()










