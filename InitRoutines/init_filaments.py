#!/usr/bin/env python

"""@package docstring
File: init_filaments.py
Author: Adam Lamson
Email: adam.lamson@colorado.edu
Description:
"""
from mpl_toolkits import mplot3d  # for vizualization
import sys
import yaml
import numpy as np
import matplotlib.pyplot as plt
import math
import argparse
from numba import jit
# %matplotlib notebook


def parseArgs():
    parser = argparse.ArgumentParser(prog='init_filaments.py')
    parser.add_argument("input", default=None,
                        help="File with initialization parameters.")
    parser.add_argument(
        '-T',
        '--tactoid',
        action='store_true',
        default=False,
        help='Generate nematic tactoid initial conditions.')
    parser.add_argument(
        '-I',
        '--isotropic',
        action='store_true',
        default=False,
        help='Generate isotropic initial conditions. Overlaps may occur.')
    opts = parser.parse_args()
    return opts


class Filament():
    def __init__(self, center, director, leng, radius, gid):
        self.center = center
        self.director = director / np.linalg.norm(director)
        self.length = leng
        self.radius = radius
        self.pos_start = center - leng / 2 * self.director
        self.pos_end = center + leng / 2 * self.director
        self.gid = gid  # defined outisde the class

    def Plot3(self, ax):
        ax.plot3D([self.pos_start[0], self.pos_end[0]], [
                  self.pos_start[1], self.pos_end[1]], [self.pos_start[2], self.pos_end[2]], 'red')

    def GetStringtoWrite(self):
        return 'C {0} {1} {2:0.6f} {3:0.6f} {4:0.6f} {5:0.6f} {6:0.6f} {7:0.6f}\n'.format(
            self.gid, self.radius,
            self.pos_start[0], self.pos_start[1], self.pos_start[2],
            self.pos_end[0], self.pos_end[1], self.pos_end[2])

    def __repr__(self):
        return "Filament()"

    def __str__(self):
        return 'Filament {0}:\n  pos_start: {1}\n  pos_end: {2}\n  length: {3}\n  radius: {4}'.format(
            self.gid, self.pos_start, self.pos_end, self.length, self.radius)


@jit
def minDistBetweenTwoFil(p1, p2, p3, p4):
    """ Adapted from matlab
     https://www.mathworks.com/matlabcentral/fileexchange/32487-shortest-distance-between-two-line-segments
     which adapted this from Dan Sunday's Geometry Algorithms originally
     written in C++
     http://softsurfer.com/Archive/algorithm_0106/algorithm_0106.htm#dist3D_Segment_to_Segment

     Computes the minimum distance between two line segments. Code
     is adapted for Matlab from Dan Sunday's Geometry Algorithms originally
     written in C++
     http://softsurfer.com/Archive/algorithm_0106/algorithm_0106.htm#dist3D_Segment_to_Segment

     Usage: Input the start and end x,y,z coordinates for two line segments.
     p1, p2 are [x,y,z] coordinates of first line segment and p3,p4 are for
     second line segment.
     Output: scalar minimum distance between the two segments.
      Example:
        P1 = [0 0 0];     P2 = [1 0 0];
        P3 = [0 1 0];     P4 = [1 1 0];
    dist = minDistBetweenTwoFil(P1, P2, P3, P4)"""

    u = p1 - p2
    v = p3 - p4
    w = p2 - p4

    a = np.dot(u, u)
    b = np.dot(u, v)
    c = np.dot(v, v)
    d = np.dot(u, w)
    e = np.dot(v, w)
    D = a * c - b * b
    sD = D
    tD = D

    SMALL_NUM = 0.00000001

    # compute the line parameters of the two closest points
    if D < SMALL_NUM:  # the lines are almost parallel
        sN = 0.0     # force using point P0 on segment S1
        sD = 1.0     # to prevent possible division by 0.0 later
        tN = e
        tD = c
    else:             # get the closest points on the infinite lines
        sN = (b * e - c * d)
        tN = (a * e - b * d)
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
    if np.absolute(sN) < SMALL_NUM:
        sc = 0.0
    else:
        sc = sN / sD

    if np.absolute(tN) < SMALL_NUM:
        tc = 0.0
    else:
        tc = tN / tD

    # get the difference of the two closest points
    dP = w + (sc * u) - (tc * v)  # = S1(sc) - S2(tc)
    distance = np.linalg.norm(dP)
    # outV = dP

    # outV = outV      # vector connecting the closest points
    # cp_1 = p2+sc*u  # Closest point on object 1
    # cp_2 = p4+tc*v  # Closest point on object 2

    return distance


def check_fil_overlap(fil, fil_list, min_dist_bound):
    min_dist = 1e6
    for fil_check in fil_list:
        check_dist = minDistBetweenTwoFil(fil.pos_start, fil.pos_end,
                                          fil_check.pos_start,
                                          fil_check.pos_end)
        min_dist = check_dist if check_dist < min_dist else min_dist
        if min_dist < min_dist_bound:
            return False
    return True


def getLengthRandomExp(rng, Lmean, Lmin=0, Lmax=1e6):
    # Define a random length sampler
    if Lmin >= Lmax:
        raise ValueError("Lmax must be greater than Lmin.")
    val = -1
    i = 0
    while val < Lmin or val > Lmax:
        val = rng.exponential(Lmean)
        i += 1
        if i > 5000:
            raise RuntimeError("Taking too long to find random number. "
                               "Try increasing difference between Lmin and Lmax.")
    return val


def sph2cart(r, phi, theta):
    # polar to cartesian coordinates
    return np.array([
        r * math.sin(theta) * math.cos(phi),
        r * math.sin(theta) * math.sin(phi),
        r * math.cos(theta)
    ])


def cart2sph(x, y, z):
    # cartesian to polar coordinates
    XsqPlusYsq = x**2 + y**2
    r = math.sqrt(XsqPlusYsq + z**2)               # r
    elev = math.atan2(z, math.sqrt(XsqPlusYsq))     # theta
    az = math.atan2(y, x)                           # phi
    return np.array([r, az, elev])


def getFilamentLengthInsideSphere(f1, rad):
    # volume of filament inside sphere

    # Find coordinates along the filament whose r coord in (r,phi, theta) lies inside the sphere.
    # Segment the filament into 100 points
    lens = np.linspace(0, f1.length, 100)
    num_in = 0
    for el in lens:
        [x, y, z] = f1.pos_start + el * f1.director
        pt = cart2sph(x, y, z)
        if pt[0] < rad:
            num_in += 1

    return (num_in / 100) * f1.length


def getRand3PointInSphere(rng, rad):
    # Get random 3d cartesian point() inside sphere
    d = 2 * rad
    while d > rad:
        xyz = rng.uniform(-rad, rad, 3)
        d = np.linalg.norm(xyz)
    return [xyz[0], xyz[1], xyz[2]]


def attemptAddFilament(center, L, director, D, f_list, gid):
    # Initialize a  filament of length L
    fil = Filament(center, director, L, D * .5, gid)

    # find minimum distance to other filaments, and accept if distance is more
    # than filament diameter
    status = check_fil_overlap(fil, f_list, D)
    if status:
        f_list.append(fil)
        gid += 1
        print(fil)
    return f_list, gid, status


def generate_nematic_sphere(fil_diam=.007, Lmean=0.18, Lmin=.09, Lmax=.028,
                            director=[0., 0., 1.],
                            sphere_radius=.64,
                            n_fil=2000,
                            use_pack_frac=False,
                            use_exp_filament_length=True,
                            **kwargs):
    """TODO: Docstring for generate_nematic_sphere.
    @return: TODO

    """
    gid = 0
    f_list = []  # list to store filaments
    vol_sphere = (4. / 3.) * math.pi * (sphere_radius**3)
    rng = np.random.default_rng()  # initialize generator instance
    if use_pack_frac:

        v_enclosed = 0  # volume occupied by filaments inside the sphere
        cpf = 0  # current packing fraction tracker
        while cpf < kwargs['pack_frac']:
            # Sample a random (x,y,z inside a sphere)
            center = getRand3PointInSphere(rng, sphere_radius)

            if use_exp_filament_length:
                L = getLengthRandomExp(rng, Lmean, Lmin, Lmax)
            else:
                L = Lmean

            f_list, gid, status = attemptAddFilament(
                center, L, director, fil_diam, f_list, gid)
            if status:
                # length of filaments inside the sphere
                len_inside = getFilamentLengthInsideSphere(
                    f_list[-1], sphere_radius)
                # volume of filament inside the sphere
                v_enclosed += math.pi * (fil_diam * .5)**2 * len_inside
                cpf = v_enclosed / vol_sphere  # current packing fraction
                print('CPF = {}'.format(cpf))
    else:

        # sample n_fil length values and sort them
        if use_exp_filament_length:
            lens = sorted([getLengthRandomExp(rng, Lmean, Lmin, Lmax)
                           for i in range(n_fil)], reverse=True)
        else:
            lens = [Lmean] * n_fil

        for idx, el in enumerate(lens):
            status = False
            failCount = 0
            while not status:
                failCount += 1
                # Sample a random (x,y,z inside a sphere)
                center = getRand3PointInSphere(rng, sphere_radius)
                f_list, gid, status = attemptAddFilament(
                    center, el, director, fil_diam, f_list, gid)
                if failCount > 5000:
                    raise Exception(
                        'Failed to add a filament with 5000 trials.\nStopping...\n...\n...\n Total Added Filaments = {}'.format(
                            len(f_list)))
    return f_list


def generate_isotropic_cube(fil_diam=.007, Lmean=0.18, Lmin=.09, Lmax=.028,
                            side_length=.64,
                            n_fil=2000,
                            use_pack_frac=False,
                            use_exp_filament_length=True,
                            **kwargs):
    """TODO: Docstring for generate_isotropic_cube.

    @param fil_diam TODO
    @param Lmean TODO
    @return: TODO

    """
    rng = np.random.default_rng()  # initialize generator instance
    f_list = []  # list to store filaments
    if use_exp_filament_length:
        lens = sorted([getLengthRandomExp(rng, Lmean, Lmin, Lmax)
                       for i in range(n_fil)], reverse=True)
    else:
        lens = [Lmean] * n_fil

    for idx, el in enumerate(lens):
        center = rng.uniform(0, side_length, 3)
        director = rng.uniform(0, 1, 3)
        director /= np.linalg.norm(director)
        f_list += [Filament(center, director, el, .5 * fil_diam, idx)]

    return f_list


def main():
    """TODO: Docstring for main.

    @param yml_file TODO
    @return: TODO

    """
    opts = parseArgs()
    with open(opts.input, 'r') as yf:
        p_dict = yaml.safe_load(yf)

    if opts.tactoid:
        print(":::Making nematic tactiod initial conditions:::")
        fil_list = generate_nematic_sphere(**p_dict)

    elif opts.isotropic:
        print(":::Making isotropic initial conditions:::")
        fil_list = generate_isotropic_cube(**p_dict)
    else:
        print("Must select type of initial conditions to generate.")
        return

    # first offset system
    fname = p_dict.get('out_file_name', 'TubuleInitial.dat')
    with open(fname, 'w') as filer:
        filer.write('# Initial configuration of rods\n#\n')
        system_offset = np.asarray(p_dict.get('system_offset', [0, 0, 0]))
        for fil in fil_list:
            fil.pos_start += np.array(system_offset)
            fil.pos_end += np.array(system_offset)
            fil.center += np.array(system_offset)
            filer.write(fil.GetStringtoWrite())


##########################################
if __name__ == "__main__":
    main()
