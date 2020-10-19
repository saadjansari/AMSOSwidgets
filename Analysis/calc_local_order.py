import numpy as np
from numba import njit
from common_func import *

@njit
def calc_local_polar_order( c, orients, boxsize):
    # Calculate the local polar order

    porder_par = 0
    porder_apar = 0

    # coordinate indices within the distance range
    q = pair_partition_func_centers(c)

    for idx in np.arange( q.shape[0]):

        o_parallel = np.zeros(3)
        o_aparallel = np.zeros(3)
        n_parallel = 0
        n_aparallel = 0

        # orientation of reference filament
        o1 = orients[idx,:]
        dp = np.dot(orients, o1)
        scaled_orients = np.multiply( orients.transpose(), q[idx,:]).transpose()
        for idx1 in np.arange(q.shape[1]):
            if dp[idx1] >0:
                o_parallel += scaled_orients[idx1,:]
                n_parallel += 1
            else:
                o_aparallel += scaled_orients[idx1,:]
                n_aparallel += 1

        if n_parallel > 0:
            porder_par += np.linalg.norm( o_parallel/n_parallel)
        if n_aparallel > 0:
            porder_apar += np.linalg.norm( o_aparallel/n_aparallel)
        pdb.set_trace()

    return porder_par/c.shape[0], porder_apar/c.shape[0]
