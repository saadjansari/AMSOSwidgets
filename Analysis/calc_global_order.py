import numpy as np
from numba import njit

def calc_nematic_order(orient_array):
    # Calculates the maximum eigenvalue of the nematic tensor Q

    if not np.any( orient_array):
        return np.nan

    # calculate Q tensor
    Q = calc_nematic_tensor(orient_array)
    S = np.sqrt(np.tensordot(Q, Q)*1.5)
    return S

# @njit
def calc_polar_order(orient_array):
    # Calculates the polar order given the array of filament orientations

    if not np.any( orient_array):
        return np.nan

    # Initialize P vector
    Pvec = np.zeros(3)

    # Take a mean of all orientations
    for irow in np.arange( orient_array.shape[0]):
        Pvec += orient_array[irow,:]
    Pvec = Pvec / orient_array.shape[0]
    P = np.linalg.norm( Pvec)
    return P
    

@njit
def calc_nematic_tensor( orient_array):
    # Calculates the nematic tensor Q

    # initialize Q tensor
    Q = np.zeros((3,3))
    # sum over all filaments, taking their outer products
    for irow in np.arange(orient_array.shape[0]):
        Q += np.outer( orient_array[irow,:], orient_array[irow,:])

    # complete mean calculation by dividing by number of filaments and subtracting identity.
    Q = Q/orient_array.shape[0] - np.identity(3)/3
    return Q

@njit
def calc_z_ordering( orient_array):
    # Calculate alignment in z

    if not np.any( orient_array):
        return np.nan

    z_axis = np.array( [0.,0.,1.])
    sum_all = 0
    for idx in np.arange( orient_array.shape[0]):
        sum_all += np.absolute( orient_array[idx,:].dot(z_axis) )
    return sum_all/orient_array.shape[0]
    
