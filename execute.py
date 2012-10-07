import numpy
from numpy.random import normal
from numpy.random import randint

import src.multics as multics


def sample_nuisance(s, sp, ss, ssp, dlog):
    return s, sp, ss, ssp, dlog


def sample_assignment(


def calc_energy():
    return 0


def calc_log_bias():
    return 0


def evaluate():

    s, sp, ss, ssp, dlog = update_variables(s, sp, ss, ssp, dlog)

    energy  = calc_energy()
    bias    = calc_log_bias()

    return energy + bias



s    = [0.4, 2.9]
sp   = [0.4, 2.9]

ss   = [1.0, 1.0]
ssp  = [1.0, 1.0]

dlog = [1.0, 1.0]


n_steps = int(1E6)

