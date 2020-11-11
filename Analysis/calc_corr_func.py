import numpy as np
from numba import njit
from coord_transform import *


def calc_cyl_corr_single( c_ref, c_all, idx_exclude, boxsize, axis_cyl=np.array([0,0,1])):
    # calculates a binned (x^2+y^2) vs z binned matrix for distributions of points relative to reference point.

    c_tall = unfold_coordinates(c_all[ np.arange(c_all.shape[0])!=idx_exclude], c_ref,boxsize)-c_ref
    # transform all other coordinates
    c_tall = transform_coords(
        c_tall,
        np.array([0,0,0]),
        axis_cyl)
    x2y2sqrt = np.sqrt( c_tall[:,0]**2 + c_tall[:,1]**2 )
    z = c_tall[:,2]
    return x2y2sqrt, z


def calc_cyl_corr( c_all, orients, boxsize, **kwargs):
    # calculates a binned (x^2+y^2) vs z binned matrix for all pairs of distributions of points


    sz = c_all.shape[0]
    xyall = np.zeros( sz*(sz-1))
    zall = np.zeros( sz*(sz-1))
    for idx in np.arange(sz):
        if 'axis_cyl' in kwargs.keys():
            new_axis = kwargs['axis_cyl']
        else:
            new_axis = orients[idx,:] 
        xyall[idx*(sz-1):idx*(sz-1)+(sz-1)], zall[idx*(sz-1):idx*(sz-1)+(sz-1)] = calc_cyl_corr_single(
                c_all[idx,:],
                c_all, idx, boxsize, new_axis)
    return xyall, zall


# radial distribution function
def calc_rdf_single( c_ref, c_all, boxsize):
    # calculates a binned (x^2+y^2) vs z binned matrix for distributions of points relative to reference point.
    c_tall = unfold_coordinates(c_all, c_ref, boxsize)-c_ref
    r = np.sqrt( c_tall[:,0]**2 + c_tall[:,1]**2 + c_tall[:,2]**2 )
    return r


def calc_rdf( c_all, boxsize):
    # calculates distances for all pairs of distributions of points

    sz = c_all.shape[0]
    rall = np.zeros( sz*(sz-1))
    for idx in np.arange(sz):
        c_aall = c_all[ np.arange(c_all.shape[0])!=idx]
        rall[idx*(sz-1):idx*(sz-1)+(sz-1)] = calc_rdf_single(
                c_all[idx,:],c_aall, boxsize)
    return rall

