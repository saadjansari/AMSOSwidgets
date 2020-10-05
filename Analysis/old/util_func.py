import numpy as np
from numba import jit
from calc_global_order import *

# Filament class
class Filament():
    def __init__(self, pos0, pos1, radius,gid):
        self.radius = radius
        self.pos0 = pos0
        self.pos1 = pos1
        self.gid = gid
    def GetCenter(self,boxsize):
        return getMean(self.pos0, self.pos1,boxsize)
    def GetLength(self,boxsize):
        xi = getDistance(self.pos1,self.pos0,boxsize)
        return np.linalg.norm( xi)
    def GetOrientation(self):
        xi = self.pos1 - self.pos0
        return xi/np.sqrt(xi.dot(xi))

    def Plot3(self,ax, col="red"):
        ax.plot3D( [self.pos0[0], self.pos1[0]], [self.pos0[1], self.pos1[1]], [self.pos0[2], self.pos1[2]], col)
    def GetStringtoWrite(self):
        return 'C {0} {1} {2:0.6f} {3:0.6f} {4:0.6f} {5:0.6f} {6:0.6f} {7:0.6f}\n'.format(
        self.gid, self.radius,
        self.pos0[0], self.pos0[1], self.pos0[2],
        self.pos1[0], self.pos1[1], self.pos1[2])
    def __repr__(self):
        return "Filament()"
    def __str__(self):
        return 'Filament {0}:\n  pos0: {1}\n  pos1: {2}\n  radius: {3}'.format(self.gid, self.pos0, self.pos1,self.radius)

# Protein class
class Protein():
    def __init__(self, pos0, pos1, link0, link1, gid):
        self.pos0 = pos0
        self.pos1 = pos1
        self.link0 = link0
        self.link1 = link1
        self.gid = gid
    def GetCenter(self,boxsize):
        return getMean(self.pos0, self.pos1,boxsize)
    def GetLength(self,boxsize):
        xi = getDistance(self.pos1,self.pos0,boxsize)
        return np.linalg.norm( xi)
    def GetOrientation(self):
        if link0 != -1 and link1 != -1:
            xi = self.pos1 - self.pos0
            return xi/np.sqrt( xi.dot(xi))
        else:
            return None

    def Plot3(self,ax,col="blue"):
        ax.plot3D( [self.pos0[0], self.pos1[0]], [self.pos0[1], self.pos1[1]], [self.pos0[2], self.pos1[2]], col)
    def GetStringtoWrite(self):
        return 'P {0} 0 {2:0.6f} {3:0.6f} {4:0.6f} {5:0.6f} {6:0.6f} {7:0.6f} {8} {9} \n'.format(
        self.gid, self.radius,
        self.pos0[0], self.pos0[1], self.pos0[2],
        self.pos1[0], self.pos1[1], self.pos1[2],
        self.link0, self.link1)
    def __repr__(self):
        return "Protein()"
    def __str__(self):
        return 'Protein {0}:\n  pos0: {1}\n  pos1: {2}\n  Links: {3}--{4}'.format(self.gid, self.pos0, self.pos1, self.link0, self.link1)

def getDistance(p0,p1,boxsize):
    # distance between two points in the nearest image convention
    # can use multidimensional arrays for distances between multiple points
    dist = np.absolute( p1-p0)
    for idx in range(dist.shape[-1]):
        if len(dist.shape) == 1:
            k = np.floor( dist[idx]/(0.5*boxsize[idx]))
            dist[idx] -= k*boxsize[idx]
        elif len(dist.shape) == 2:
            k = np.floor( dist[:,idx]/(0.5*boxsize[idx]))
            dist[:,idx] -= k*boxsize[idx]
        elif len(dist.shape) == 3:
            k = np.floor( dist[:,:,idx]/(0.5*boxsize[idx]))
            dist[:,:,idx] -= k*boxsize[idx]
    return np.absolute(dist)

def getMean(p0,p1,boxsize):
    # mean of the two points in the nearest image
    dist = np.absolute(p1-p0)
    for idx in range(dist.shape[-1]):
        if len(dist.shape) == 1:
            k = np.floor( dist[idx]/(0.5*boxsize[idx]))
            p1[idx] -= k*boxsize[idx]
        elif len(dist.shape) == 2:
            k = np.floor( dist[:,idx]/(0.5*boxsize[idx]))
            p1[:,idx] -= k*boxsize[idx]
        elif len(dist.shape) == 3:
            k = np.floor( dist[:,:,idx]/(0.5*boxsize[idx]))
            p1[:,:,idx] -= k*boxsize[idx]
    return (p0 + p1)/2


def getCIBootstrap(data,alpha,n_iter=10000):
    boot_mu = [np.mean(np.random.choice(x, len(x))) for _ in range(iteration)]
    beta = (1-alpha)/2
    lo_x_boot = np.percentile(boot_mu, 100*beta)
    hi_x_boot = np.percentile(boot_mu, 100*(1-beta))
    return (lo_x_boot, hi_x_boot)

class AnalysisBook():
    def __init__(self,boxsize):
        self.boxsize = boxsize
    # refpair to all other filament distances
    def R2Adistance(self, rfil, fils):
        dist = np.zeros(( len(fils)-1,3))
        cnt=0
        for cnt,f1 in enumerate(fils):
            if f1.gid is not rfil.gid:
                dist[cnt,:] = getDistance( rfil.GetCenter(self.boxsize), f1.GetCenter(self.boxsize),self.boxsize)
                cnt+=1
        return dist
    # pair to pair distances
    def P2Pdistance(self, fils):
        # only computer for upper triangle. Then add the transpose to get bottom triangle
        p00 = np.zeros(( len(fils), len(fils),3))
        p01 = np.zeros(( len(fils), len(fils),3))
        p10 = np.zeros(( len(fils), len(fils),3))
        p11 = np.zeros(( len(fils), len(fils),3))
        for i1,f1 in enumerate(fils):
            for i2,f2 in enumerate(fils):
                if i1 > i2:
                    p00[i1,i2,:] = f1.pos0
                    p01[i1,i2,:] = f1.pos1
                    p10[i1,i2,:] = f2.pos0
                    p11[i1,i2,:] = f2.pos0
        c0 = getMean( p00,p01, self.boxsize)
        c1 = getMean( p10,p11, self.boxsize)
        dist = getDistance(c0,c1,self.boxsize)
        dist = dist + np.transpose(dist, (1,0,2) )
        return dist
    # get center of mass
    def COM(self, fils):
        sumList = np.zeros((1,3))
        for fil in fils:
            sumList += fil.GetLength(self.boxsize)*getDistance(np.array([0,0,0]), fil.GetCenter(self.boxsize), self.boxsize )
        COM = sumList / np.sum( [f.GetLength(self.boxsize) for f in fils])
        return COM[0]
    # Get reference frame nematic basis vectors
    def GetBasisVectors( self,fils):
        # basis vectors are the eigenvectors of the nematic tensor Q in order of magnitude of their eigenvalues
        Q,_ = self.NematicTensor(fils)
        w,v = np.linalg.eig( Q)
        idx = np.argsort(np.sqrt(w*w))
        v = v[:,idx]
        w = w[idx]
        idx = [0,1,2]
        if np.any( np.around( np.cross( v[:,0], v[:,1]),2) != np.around( v[:,2],2) ):
            idx1 = idx[0]
            idx2 = idx[1]
            idx = [idx2, idx1, idx[2]]
        w = w[idx]
        v = v[:,idx]
        return w,v
    def FilsIdxInsideTac_T1_T2( self, tilist,t1,t2):
    # find the filaments that are remain inside a tactoid between times t1 and t2
    # the ids of filaments inside tac are stored in tilist
        cc = tilist[t1]
        for tt in range(t1, t2+1):
            cc = list( set.intersection(set(cc), set(tilist[tt])))
            #cc = [val for val in tilist[tt] if val in cc]
        return cc
    def FilsIdxOutsideTac_T1_T2( self, tilist,t1,t2, nF):
    # find the filaments that remain outside a tactoid between times t1 and t2
    # the ids of filaments inside tac are stored in tilist
        cc = list(range(0,nF))
        for tt in range(t1, t2+1):
            for xx in tilist[tt]:
                if xx in cc:
                    cc.remove(xx)
        return cc
    def FilsInsideTac( self, flist, idx):
    # find the filaments that are remain inside a tactoid between times t1 and t2
    # the ids of filaments inside tac are stored in tilist
        ff = []
        for ec in idx:
            ff.append( flist[ec])
        return ff
    def P2PInsideTac( self, p2p, idx):
        # returns a mini matrix corresponding to the pairs inside tac. Needs the full p2p matrix.

        # delete rows and column not corresponding to those in idx
        rm = list( set.difference(set( range( p2p.shape[0]) ), set( idx )))
        #rm = [xx for xx in range( p2p.shape[0]) if xx not in idx]
        p2pdel = np.delete( p2p, rm, axis=0)
        p2pdel = np.delete( p2pdel, rm, axis=1)
        return p2pdel

def PairSeparation( p2pmat):
    pp = np.sum( p2pmat**2, axis=2)
    pp = np.sqrt( pp[ np.triu_indices(pp.shape[0], k = 1)] )
    mu = np.mean(pp)
    stderr = np.std(pp)/np.sqrt(pp.size-1)
    return mu,mu-2*stderr,mu+2*stderr


@jit( nopython=True)
def dist2all(c,nT,nF,boxsize,boolRodsInside):
    # Calculate mean pair-pair filament distances for all filaments
    ppt=np.zeros(nF)
    p2p = np.zeros((nF,nF))

    for jf in range(nF):
        if not boolRodsInside[jf]:
            continue
        res = 0
        dist = np.absolute( c[jf,:]-c)
        # get distance to all other filaments
        for jf2 in range(nF):
            if not boolRodsInside[jf2] or jf==jf2:
                continue
            elif p2p[jf2,jf] == 0:
                res2=0
                # get distance
                #dist = np.absolute( c[jf,:]-c[jf2,:])
                for idx in range(3):
                    k = np.floor( dist[jf2,idx]/(0.5*boxsize[idx]))
                    dist[jf2,idx] -= k*boxsize[idx]
                    res2 += dist[jf2,idx]**2
                p2p[jf,jf2] = np.sqrt(res2)
                res += p2p[jf,jf2]
            else:
                res += p2p[jf2,jf]
        ppt[jf] = res/(np.sum(boolRodsInside)-1)
    return ppt

@jit( nopython=True)
def meanPairSep(c,nT,nF,boxsize,boolRodsInside):
    # Calculate and store mean pair-pair filament distances for all filaments for all time
    mu = np.zeros(nT)
    sig = np.zeros(nT)
    for it in range(nT):
        ppt=dist2all(c[it,:,:],nT,nF,boxsize,boolRodsInside[it])
        mu[it] = np.mean(ppt[ boolRodsInside[it] ])
        sig[it] = np.std(ppt[ boolRodsInside[it] ])
    return mu,sig

@jit( nopython=True)
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

# @jit( nopython=True)
# def within_range(c,nF,dist_max,boxsize,boolRodsInside):
#     # Calculate mean pair-pair filament distances for all filaments
#     p2p = np.zeros((nF,nF))
#     for jf in range(nF):
#         if not boolRodsInside[jf]:
#             continue
#         dist = np.absolute( c[jf,:]-c)
#         # get distance to all other filaments
#         for jf2 in range(nF):
#             if not boolRodsInside[jf2]:
#                 continue
#             elif jf==jf2:
#                p2p[jf,jf2] = 100*dist_max
#                p2p[jf2,jf] = 100*dist_max
#             elif p2p[jf2,jf] == 0:
#                 res2=0
#                 # get distance
#                 for idx in range(3):
#                     k = np.floor( dist[jf2,idx]/(0.5*boxsize[idx]))
#                     dist[jf2,idx] -= k*boxsize[idx]
#                     res2 += dist[jf2,idx]**2
#                 p2p[jf,jf2] = np.sqrt(res2)
#             else:
#                 p2p[jf,jf2] =p2p[jf2,jf]
#     return (p2p < dist_max)
