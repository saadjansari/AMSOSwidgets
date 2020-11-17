import yaml
import numpy as np
import networkx as nx
from pathlib import Path
from numba import njit,jit
import matplotlib.pyplot as plt

# Specify filenames
parent_path = Path('.')
cfg_old = parent_path / 'RunConfig_old.yaml'
cfg_new = parent_path / 'RunConfig.yaml'
tub_old = parent_path / 'TubuleOld.dat'
prot_old = parent_path / 'ProteinOld.dat'
tub_new = parent_path / 'TubuleInitial.dat'
prot_new = parent_path / 'ProteinInitial.dat'

write_file = True
run_visualization = False

# Extract simbox size from yaml
def simbox_size_from_cfg(cfg_path):
    with open(cfg_path) as ff:
        dat = yaml.safe_load(ff)

    return np.array(dat['simBoxLow']),np.array(dat['simBoxHigh'])

def calc_distance_pbc(p0,p1,boxsize):
    # distance between two points in the nearest image convention
    # can use multidimensional arrays for distances between multiple points
    dist = np.absolute( p1-p0)
    for idx in np.arange(dist.shape[-1]):
        if len(dist.shape) == 1:
            k = np.floor( dist[idx]/(0.51*boxsize[idx]))
            dist[idx] -= k*boxsize[idx]
        elif len(dist.shape) == 2:
            k = np.floor( dist[:,idx]/(0.51*boxsize[idx]))
            dist[:,idx] -= k*boxsize[idx]
        # elif len(dist.shape) == 3:
            # k = np.floor( dist[:,:,idx]/(0.5*boxsize[idx]))
            # dist[:,:,idx] -= k*boxsize[idx]
    return np.absolute(dist)

def get_rand_point(low,high):
    pt = []
    for lb,ub in zip(low,high):
        pt.append( np.random.uniform(lb,ub) )
    return np.array( pt)

# Filament class
class Filament():
    def __init__(self, pos0, pos1, radius,gid):
        self.radius = radius
        self.pos0 = pos0
        self.pos1 = pos1
        self.gid = gid
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
    def GetOrientation(self):
        if link0 != -1 and link1 != -1:
            xi = self.pos1 - self.pos0
            return xi/np.sqrt( xi.dot(xi))
        else:
            return None

    def Plot3(self,ax,col="blue"):
        ax.plot3D( [self.pos0[0], self.pos1[0]], [self.pos0[1], self.pos1[1]], [self.pos0[2], self.pos1[2]], col)
    def GetStringtoWrite(self):
        return 'P {0} 0 {1:0.6f} {2:0.6f} {3:0.6f} {4:0.6f} {5:0.6f} {6:0.6f} {7} {8} \n'.format(
        self.gid,
        self.pos0[0], self.pos0[1], self.pos0[2],
        self.pos1[0], self.pos1[1], self.pos1[2],
        self.link0, self.link1)
    def __repr__(self):
        return "Protein()"
    def __str__(self):
        return 'Protein {0}:\n  pos0: {1}\n  pos1: {2}\n  Links: {3}--{4}'.format(self.gid, self.pos0, self.pos1, self.link0, self.link1)

def check_fil_overlap(fil, fil_list, min_dist_bound):
    min_dist = 1e6
    for fil_check in fil_list:
        check_dist = minDistBetweenTwoFil(fil.pos0, fil.pos1,
                                          fil_check.pos0, fil_check.pos1)
        min_dist = check_dist if check_dist < min_dist else min_dist
        if min_dist < min_dist_bound:
            return False
    return True

# Find Minimum Distance between two filaments
@jit
def minDistBetweenTwoFil(p1, p2, p3, p4):
    # Adapted from matlab
    # https://www.mathworks.com/matlabcentral/fileexchange/32487-shortest-distance-between-two-line-segments
    # which adapted this from Dan Sunday's Geometry Algorithms originally written in C++
    # http://softsurfer.com/Archive/algorithm_0106/algorithm_0106.htm#dist3D_Segment_to_Segment

    # p1 = f1.pos_start
    # p2 = f1.pos_end
    # p3 = f2.pos_start
    # p4 = f2.pos_end
    # Computes the minimum distance between two line segments. Code
    # is adapted for Matlab from Dan Sunday's Geometry Algorithms originally
    # written in C++
    # http://softsurfer.com/Archive/algorithm_0106/algorithm_0106.htm#dist3D_Segment_to_Segment
    # Usage: Input the start and end x,y,z coordinates for two line segments.
    # p1, p2 are [x,y,z] coordinates of first line segment and p3,p4 are for
    # second line segment.
    # Output: scalar minimum distance between the two segments.
    #  Example:
    #	P1 = [0 0 0];     P2 = [1 0 0];
    #   P3 = [0 1 0];     P4 = [1 1 0];
    #	dist = minDistBetweenTwoFil(P1, P2, P3, P4)
    #p1 = f1.pos_start
    #p2 = f1.pos_end
    #p3 = f2.pos_start
    #p4 = f2.pos_end

    u = p1 - p2
    v = p3 - p4
    w = p2 - p4

    a = np.dot(u,u)
    b = np.dot(u,v)
    c = np.dot(v,v)
    d = np.dot(u,w)
    e = np.dot(v,w)
    D = a*c - b*b
    sD = D
    tD = D

    SMALL_NUM = 0.00000001

    # compute the line parameters of the two closest points
    if D < SMALL_NUM: # the lines are almost parallel
        sN = 0.0     # force using point P0 on segment S1
        sD = 1.0     # to prevent possible division by 0.0 later
        tN = e
        tD = c
    else:             # get the closest points on the infinite lines
        sN = (b*e - c*d)
        tN = (a*e - b*d)
        if sN < 0.0:   # sc < 0 => the s=0 edge is visible
            sN = 0.0
            tN = e
            tD = c
        elif sN > sD:  # sc > 1 => the s=1 edge is visible
            sN = sD
            tN = e + b
            tD = c

    if tN < 0.0:            # tc < 0 => the t=0 edge is visible
        tN = 0.0
        # recompute sc for this edge
        if -d < 0.0:
            sN = 0.0
        elif -d > a:
            sN = sD
        else:
            sN = -d
            sD = a
    elif tN > tD:       # tc > 1 => the t=1 edge is visible
        tN = tD
        # recompute sc for this edge
        if -d + b < 0.0:
            sN = 0
        elif -d + b > a:
            sN = sD
        else:
            sN = -d + b
            sD = a

    # finally do the division to get sc and tc
    if  np.absolute(sN) < SMALL_NUM:
        sc = 0.0
    else:
        sc = sN / sD

    if np.absolute(tN) < SMALL_NUM:
        tc = 0.0
    else:
        tc = tN / tD

    # get the difference of the two closest points
    dP = w + (sc * u) - (tc * v);  # = S1(sc) - S2(tc)
    distance = np.linalg.norm(dP);
    outV = dP;

    # outV = outV      # vector connecting the closest points
    # cp_1 = p2+sc*u  # Closest point on object 1
    # cp_2 = p4+tc*v  # Closest point on object 2

    return distance

# get box sizes
sb_low,sb_high = simbox_size_from_cfg(cfg_new)
sb_low_old,sb_high_old = simbox_size_from_cfg(cfg_old)

# find filaments
flist = []
with open(tub_old, 'r') as file1:
    for line in file1:
        if line.startswith('C'):
            data = line.split()
            gid = int(data[1])
            radius = float(data[2])
            pos0 = np.array([float(data[3]), float(data[4]), float(data[5])])
            pos1 = np.array([float(data[6]), float(data[7]), float(data[8])])
            flist.append( Filament(pos0, pos1, radius,gid))


plist = []
with open(prot_old, 'r') as file2:
    for line in file2:
        if line.startswith('P'):
            data = line.split()
            gid = int(data[1])
            pos0 = np.array([float(data[3]), float(data[4]), float(data[5])])
            pos1 = np.array([float(data[6]), float(data[7]), float(data[8])])
            link0 = int(data[9])
            link1 = int(data[10])
            plist.append( Protein(pos0, pos1, link0, link1, gid))

# Get filament indices that are outside sim box
foutside = []
faccept = []
fil_lens = []
for fil in flist:
    pos0_out = np.any(fil.pos0 > sb_high) or np.any(fil.pos0 < sb_low)
    pos1_out = np.any(fil.pos1 > sb_high) or np.any(fil.pos1 < sb_low)
    if pos0_out or pos1_out:
        foutside.append( fil)
        fil_lens.append( np.linalg.norm( calc_distance_pbc(fil.pos0, fil.pos1, sb_high_old) ) )
    else:
        faccept.append(fil)

# Get protein indices that are outside sim box
idx_prot_out = []
for idx, prot in enumerate(plist):
    pos0_out = np.any(prot.pos0 > sb_high) or np.any(prot.pos0 < sb_low)
    pos1_out = np.any(prot.pos1 > sb_high) or np.any(prot.pos1 < sb_low)
    if pos0_out or pos1_out:
        idx_prot_out.append( idx)

# Loop over outside filaments and resample their positions
for fil,el in zip(foutside, fil_lens):

    # sample new position
    sample_again = True
    while sample_again:
        fil.pos0 = get_rand_point(sb_low,sb_high)
        # get ort at random
        ort = np.random.randn(3)
        ort = ort/ np.linalg.norm(ort)

        # get pos1
        fil.pos1 = fil.pos0 + ort*el
        sample_again = not check_fil_overlap(fil, faccept, 2*fil.radius)
    faccept.append(fil)

# Loop over outside proteins and resample their positions
for idx in idx_prot_out:

    # sample new position
    sample_again = True
    plist[idx].pos = get_rand_point(sb_low,sb_high)
    plist[idx].link0 = -1
    plist[idx].link1 = -1

if write_file:

    # first offset system
    with open(tub_new, "w") as filer:
        filer.write('# Initial configuration of rods\n#\n')
        for fil in flist:
            filer.write( fil.GetStringtoWrite() )

    # first offset system
    with open(prot_new, "w") as filer:
        filer.write('# Initial configuration of proteins\n#\n')
        for prot in plist:
            filer.write( prot.GetStringtoWrite() )
