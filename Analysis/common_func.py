import numpy as np
from numba import njit

@njit
def calc_distance_pbc(p0,p1,boxsize):
    # distance between two points in the nearest image convention
    # can use multidimensional arrays for distances between multiple points
    dist = np.absolute( p1-p0)
    for idx in np.arange(dist.shape[-1]):
        if len(dist.shape) == 1:
            k = np.floor( dist[idx]/(0.5*boxsize[idx]))
            dist[idx] -= k*boxsize[idx]
        elif len(dist.shape) == 2:
            k = np.floor( dist[:,idx]/(0.5*boxsize[idx]))
            dist[:,idx] -= k*boxsize[idx]
        # elif len(dist.shape) == 3:
            # k = np.floor( dist[:,:,idx]/(0.5*boxsize[idx]))
            # dist[:,:,idx] -= k*boxsize[idx]
    return np.absolute(dist)

@njit
def calc_mean_pbc(p0,p1,boxsize):
    # mean of the two points in the nearest image
    dist = np.absolute(p1-p0)
    for idx in np.arange(dist.shape[-1]):
        if len(dist.shape) == 1:
            k = np.floor( dist[idx]/(0.5*boxsize[idx]))
            p1[idx] -= k*boxsize[idx]
        elif len(dist.shape) == 2:
            k = np.floor( dist[:,idx]/(0.5*boxsize[idx]))
            p1[:,idx] -= k*boxsize[idx]
        # elif len(dist.shape) == 3:
            # k = np.floor( dist[:,:,idx]/(0.5*boxsize[idx]))
            # p1[:,:,idx] -= k*boxsize[idx]
    return (p0 + p1)/2

@njit
def unfold_coordinates(crds,c_ref,boxsize):
    # unfolded crds via the nearest image convention

    # reference coordinate
    dist = crds-c_ref
    for idx in np.arange(crds.shape[-1]):
        k = np.sign( dist[:,idx]) * np.floor( np.absolute(dist[:,idx])/(0.5*boxsize[idx]))
        crds[:,idx] -= k*boxsize[idx]
    return crds


@njit
def calc_mean_separation(coords,boxsize):
    # Calculate mean pair-pair separation

    if not coords or coords.shape[0] < 2:
        return np.nan

    # Initialize
    num_coords = coords.shape[0]
    mu = 0 
    sig = 0

    mean_p2p = calc_mean_p2p_dist( coords, boxsize)
    for idx in np.arange( num_coords):
        mu += mean_p2p[idx]
    mu = mu/num_coords
    sig = np.std( mean_p2p)
    return mu,sig


@njit
def calc_mean_p2p_dist(coords, boxsize):
    # Calculate mean separation of every filament from all other pairs

    # Initialize
    num_coords = coords.shape[0]
    mean_p2p = np.zeros( num_coords)

    # Calculate p2p distance matrix
    p2p_mat = pair_partition_func_centers(coords, boxsize)
    for idx in np.arange( num_coords):
        for idx2 in np.arange(num_coords):
            mean_p2p[idx] += p2p_mat[idx,idx2]
        mean_p2p[idx] = mean_p2p[idx]/num_coords

    return mean_p2p 


@njit
def calc_p2p_dist_mat(coords, boxsize):

    # Initialize
    num_coords = coords.shape[0]
    p2p = np.zeros( (num_coords, num_coords))

    # Loop over reference filament. For ref filament, 
    for jf in np.arange( num_coords):

        dist_sum = 0
        dist = np.absolute( coords[jf,:]-coords)

        # get distance to all other filaments
        for jf2 in np.arange( num_coords):

            # Skip if same as the reference filament
            if jf==jf2:
                continue
            # If distance with j1 and j2 switched already calculated, copy from matrix
            elif p2p[jf2,jf] != 0:
                p2p[jf,jf2] = p2p[jf2,jf]
                dist_sum += p2p[jf2,jf]
            # Otherwise, calculate distance
            else:
                # We need to find the norm distance. We will code up a manual
                # summation of squares and then take a square root. This is
                # faster with njit.
                # We apply box-size considerations to find the distance to the closest image
                dist_xyz = 0
                # get distance
                for idx in np.arange( len(boxsize)):
                    k = np.floor( dist[jf2,idx]/(0.5*boxsize[idx]))
                    dist[jf2,idx] -= k*boxsize[idx]
                    dist_xyz += dist[jf2,idx]**2
                p2p[jf,jf2] = np.sqrt( dist_xyz)
                dist_sum += p2p[jf,jf2]
    return p2p


@njit
def msd_pairs(lags,c,nT,nF,boxsize,boolRodsInside):
    mu = np.zeros(lags.shape)
    sig = np.zeros(lags.shape)
    for idx in range(lags.size):
        pptdiff = np.zeros(nT-lags[idx])
        for it in range(pptdiff.size):

            #find rods inside at both times
            rods_in = np.logical_and( boolRodsInside[it],boolRodsInside[it+lags[idx]])

            ppt1=dist2all(c[it,:,:],nT,nF,boxsize,rods_in)
            ppt2=dist2all(c[it+lags[idx],:,:],nT,nF,boxsize,rods_in)

            res = 0
            for idx2 in range(ppt1.size):
                res+=(ppt2[idx2]-ppt1[idx2])**2
            pptdiff[it] = res/ppt1.size
        mu[idx] = np.mean(pptdiff)
        sig[idx] = np.std(pptdiff)
    return mu,sig

@njit
def pair_partition_func_centers(centers, boxsize):
    
    q = np.zeros((centers.shape[0], centers.shape[0]))
    # for each pair of centers, calculate pair_partition function
    for idx1 in np.arange( centers.shape[0]):
        q[idx1,:] = pair_partition_func_i_centers( 
            centers[idx1,:], centers, q[idx1,:], boxsize)
    return np.sqrt(q)

@njit
def pair_partition_func_i_centers(r0,r1,q, boxsize):
    
    # distance between centers
    dist = calc_distance_pbc(r0, r1, boxsize)
    for idx in np.arange( dist.shape[1]):
        q+= dist[:,idx]**2
    return q
