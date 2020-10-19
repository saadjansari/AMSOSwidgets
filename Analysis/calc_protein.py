import numpy as np
from numba import njit
from common_func import *

@njit
def calc_protein_energy( lengths, rest_length):

    if not lengths:
        return np.nan, np.nan

    # Energy
    sum_energy = 0
    for ix in np.arange( len(lengths) ):
        sum_energy += 0.5*( lengths[ix]-rest_length)**2 
    xlink_mu = sum_energy / lengths.shape[0]

    sum_energy_sq = 0
    for ix in np.arange( len(lengths) ):
        sum_energy_sq += ( 0.5*( lengths[ix]-rest_length)**2 - xlink_mu)
    xlink_std = np.sqrt( sum_energy_sq/lengths.shape[0])

    return xlink_mu, xlink_std
