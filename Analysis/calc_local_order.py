import numpy as np
from numba import njit

def calcLocalPolarOrder( c, boxsize):

    dist_max = 4*0.025

    dist_within = within_range(c,nF,dist_max, boxsize)
    porder_par = []
    porder_apar = []
    for ifil in range(nF):
        # fils_neatest = ftlist[it][ np.argwhere( dist_within[ifil,:]).transpose()[0].tolist()]
        fils_nearest = list( compress( , dist_within[ifil,:].tolist()))
        if not fils_nearest:
            continue
            orient_ref = ftlist[it][ifil].GetOrientation( boxsize)
            orientList_par = []
            orientList_apar = []
            for fil in fils_nearest:
                orient_new = fil.GetOrientation( boxsize)
                if orient_ref.dot( orient_new) / (1*1) >= 0:
                    orientList_par.append( orient_new)
                else:
                    orientList_apar.append( orient_new)
            PList_par = np.array(orientList_par)
            PList_apar = np.array(orientList_apar)
            if PList_par.any():
                porder_par.append(np.linalg.norm( np.mean(PList_par, axis=0)) )
            if PList_apar.any():
                porder_apar.append(np.linalg.norm( np.mean(PList_apar, axis=0)) )
        if porder_par:
            Porder_par[:,it] = [np.mean( porder_par), np.std(porder_par)]
        if porder_apar:
            Porder_apar[:,it] = [np.mean( porder_apar), np.std(porder_apar)]

    return Porder_par, Porder_apar, Porder_par_cluster, Porder_apar_cluster

# }}}

@njit
def within_range(c,nF,dist_max,boxsize):
    # Calculate mean pair-pair filament distances for all filaments
    p2p = np.zeros((nF,nF))
    for jf in range(nF):
        dist = np.absolute( c[jf,:]-c)
        # get distance to all other filaments
        for jf2 in range(nF):
            if jf==jf2:
               p2p[jf,jf2] = 100*dist_max
               p2p[jf2,jf] = 100*dist_max
            elif p2p[jf2,jf] == 0:
                res2=0
                # get distance
                for idx in range(3):
                    k = np.floor( dist[jf2,idx]/(0.5*boxsize[idx]))
                    dist[jf2,idx] -= k*boxsize[idx]
                    res2 += dist[jf2,idx]**2
                p2p[jf,jf2] = np.sqrt(res2)
            else:
                p2p[jf,jf2] =p2p[jf2,jf]
    return (p2p < dist_max)
