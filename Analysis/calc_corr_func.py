import numpy as np
from numba import njit
from coord_transform import *

def calc_cyl_corr_single( c_ref, c_all, axis_cyl):
    # calculates a binned (x^2+y^2) vs z binned matrix for distributions of points relative to reference point.

    # transform all other coordinates
    c_tall = transform_coords( 
            c_all-c_ref, 
            np.array([0,0,0]),
            axis_cyl)
    
    x2y2sqrt = np.sqrt( 
            c_tall[:,0]**2 + c_tall[:,1]**2 )
    z = c_tall[:,2]

    return x2y2, z

def calc_cyl_corr( c_all, axis_cyl, bins):
    # calculates a binned (x^2+y^2) vs z binned matrix for all pairs of distributions of points

    sz = c_all.shape[0]
    xyall = np.zeros( sz**2)
    zall = np.zeros( sz**2)
    for idx in np.arange(sz):
        x2y2[idx*sz:idx*sz+sz-1], zidx*sz:idx*sz+sz-1] = calc_cyl_corr_single( 
                c_all[idx,:],
                c_all, axis_cyl)






    
    

