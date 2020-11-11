#!/usr/bin/env python3

"""
Created on Fri Sep  4 18:00:20 2020

@author: Saad J. Ansari

Description: This script initializes filaments in a random orientation inside
a box of specified size. We can also enforce a cylindrical boundary oriented
in the Z direction. It can also visualize this and output a .dat file for
use as an initial TubuleInitial.dat file in the AMSOS software.

Update: Generate a ProteinInitial.dat file by simulating proteins.
"""

import numpy as np
import matplotlib.pyplot as plt
import math
# from mpl_toolkits.mplot3d import Axes3D

# =============================================================================
# Initialization -- Define constants for this project
# =============================================================================

# box_size = np.array([3.17, 3.17, 3.17])
box_size = np.array([4,4,4])

# Filament number / packing fraction
number_of_filaments = 11800
use_packing_fraction = False
packing_fraction = 0.05

# Protein number
initialize_proteins = False
number_of_proteins = 400
scale_proteins = 2

# Diameter and length of filaments (can choose to sample length from an
# exponential distribution).
diameter_of_filament = 0.007
length_mean = 0.18
sample_len_form_exp_dist = True
length_min = 0.09
length_max = 0.28

# Orientation of filaments (random or specified)
use_random_orientation = True
orientation_of_filament = np.array([0, 0, 1])

# Special Boundary
boundary = {
    'active': False,
    'type': 'Tube',
    'radius': 0.25,
    'center': np.array([2, 2, 2])
    }
# boundary = {
#     'active': True,
#     'type': 'Sphere',
#     'radius': 0.5,
#     'center': np.array([1, 1, 1])
#     }

do_file_write = True
do_visualization = False

# File name to write filament information to
fname = './TubuleInitial_cond_n6400.dat'
pname = './ProteinInitial.dat'

# Initialize a random number generator
rng = np.random.default_rng()

# =============================================================================
# Define classes and functions
# =============================================================================


# Filament class
class Filament():
    """
    A class to represent a filament.

    ...

    Attributes
    ----------
    center : vector of size (3,)
        the position of the center of the filament.
    director : vector of size (3,)
        the orientation of this filament.
    length : float
        the length of this filament.
    radius : float
        the radius of this filament.
    gid: int
        global ID
    pos_start: vector of size (3,)
        start position
    pos_end: vector of size (3,)
        end position

    Methods
    -------
    info(additional=""):
        Prints the person's name and age.
    """

    def __init__(self, center, director, length, radius, gid):
        """
        Initialize an object of class filament.

        Parameters
        ----------
        center : array of shape (3,)
            the position of the center of the filament.
        director : array of shape (3,)
            the orientation of this filament.
        length : float
            the length of this filament.
        radius : float
            the radius of this filament.
        gid : int
            global ID.

        Returns
        -------
        None.

        """
        self.center = center
        self.director = director/np.linalg.norm(director)
        self.length = length
        self.radius = radius
        self.pos_start = center - length/2 * self.director
        self.pos_end = center + length/2 * self.director
        self.gid = gid

    def plot_3D(self, ax):
        """
        Plot this filament on the 3D axes provided.

        Parameters
        ----------
        ax : axes
            3D axes to plot the filament on.

        Returns
        -------
        None.

        """
        ax.plot3D(
            [self.pos_start[0], self.pos_end[0]],
            [self.pos_start[1], self.pos_end[1]],
            [self.pos_start[2], self.pos_end[2]], 'red')

    def get_string_to_write(self):
        """
        Get a formatted string to write to a .dat file.

        Returns
        -------
        str
            Formatted string with gid, radius, start position, and end
            position.

        """
        f_string = ('C {0} {1} {2:0.6f} {3:0.6f} {4:0.6f} {5:0.6f}'
                    ' {6:0.6f} {7:0.6f}\n')

        return f_string.format(
            self.gid,
            self.radius,
            self.pos_start[0],
            self.pos_start[1],
            self.pos_start[2],
            self.pos_end[0],
            self.pos_end[1],
            self.pos_end[2])

    def __repr__(self):
        """
        Representation of this object.

        Returns
        -------
        str
            class name.

        """
        return "Filament()"

    def __str__(self):
        """
        Hooks the python print representation for this object.

        Returns
        -------
        str
            prints a description/details of this filament.

        """
        p_string = 'Filament {gid}:\n  pos_start: {pos_start}\n '
        ' pos_end: {pos_end}\n  length: {length}\n  radius: {radius}'

        return p_string.format(
            gid=self.gid,
            pos_start=self.pos_start,
            pos_end=self.pos_end,
            length=self.length,
            radius=self.radius)


# Protein class
class Protein():
    """
    A class to represent a protein.

    ...

    Attributes
    ----------
    center : vector of size (3,)
        the position of the center of the filament.
    gid : int
        Global ID.

    Methods
    -------
    info(additional=""):
        Prints the person's name and age.
    """

    def __init__(self, center, gid):
        """
        Initialize an object of class protein.

        Parameters
        ----------
        center : array of shape (3,)
            the position of the center of the filament.
        gid : int
            global ID.

        Returns
        -------
        None.

        """
        self.center = center
        self.gid = gid

    def plot_3D(self, ax):
        """
        Plot this filament on the 3D axes provided.

        Parameters
        ----------
        ax : axes
            3D axes to plot the filament on.

        Returns
        -------
        None.

        """
        ax.scatter(self.center[0], self.center[1], self.center[2], s=5, c='b')

    def get_string_to_write(self):
        """
        Get a formatted string to write to a .dat file.

        Returns
        -------
        str
            Formatted string with gid, radius, start position, and end
            position.

        """
        f_string = ('C {0} {1} {2:0.6f} {3:0.6f} {4:0.6f} {5:0.6f}'
                    '{6:0.6f} {7:0.6f}\n')

        return f_string.format(
            self.gid, 0,
            self.center[0],
            self.center[1],
            self.center[2],
            self.center[0],
            self.center[1],
            self.center[2])

    def __repr__(self):
        """
        Representation of this object.

        Returns
        -------
        str
            class name.

        """
        return "Protein()"

    def __str__(self):
        """
        Hooks the python print representation for this object.

        Returns
        -------
        str
            prints a description/details of this protein.

        """
        p_string = 'Protein {gid}:\n  center: {center}'

        return p_string.format(
            gid=self.gid,
            center=self.center)


# Minimum Distance between two filaments
def minimum_distance_between_filaments(f1, f2):
    """
    Compute mimimum distance between two filament objects.

    Adapted from matlab.
    https://www.mathworks.com/matlabcentral/fileexchange/32487-shortest-distance-between-two-line-segments
    which adapted this from Dan Sunday's Geometry Algorithms originally
    written in C++.
    http://softsurfer.com/Archive/algorithm_0106/algorithm_0106.htm#dist3D_Segment_to_Segment
    Computes the minimum distance between two line segments. Code
    is adapted for Matlab from Dan Sunday's Geometry Algorithms originally
    written in C++
    http://softsurfer.com/Archive/algorithm_0106/algorithm_0106.htm#dist3D_Segment_to_Segment
    p1 = f1.pos_start
    p2 = f1.pos_end
    p3 = f2.pos_start
    p4 = f2.pos_end

    Parameters
    ----------
    f1 : Filament object
    f2 : Filament object
        DESCRIPTION.

    Returns
    -------
    distance : float
        minimum distance between the two filaments.

    """
    p1 = f1.pos_start
    p2 = f1.pos_end
    p3 = f2.pos_start
    p4 = f2.pos_end
    u = p1 - p2
    v = p3 - p4
    w = p2 - p4
    a = np.dot(u, u)
    b = np.dot(u, v)
    c = np.dot(v, v)
    d = np.dot(u, w)
    e = np.dot(v, w)
    D = a*c - b*b
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


# Random exp length
def rand_exponential_length(rng, len_mean, len_min, len_max):
    """
    Get a random length sampled from an exponential distribution.

    Parameters
    ----------
    rng : random number generator
        rand num gen
    len_mean : float
        mean of the distribution.
    len_min : float
        lower cutoff for this distribution.
    len_max : float
        upper cutoff for this distribution.

    Returns
    -------
    val : float
        random exponential length.

    """
    val = -1  # ensure while loop is entered on first try
    while (val < len_min or
           val > len_max):
        val = rng.exponential(len_mean)
    return val


# Minimum Distance from one reference filament to all other filaments
def minimum_distance_to_all_filaments(f_ref, f_list):
    """
    Calculate minimum dist b/w a reference filament and all other filaments.

    Parameters
    ----------
    f_ref : Filament object
    f_list : list
        List of Filament objects.

    Returns
    -------
    float
        minimum distance to any other filament.

    """
    # calculate distance between this ref filament to all other filaments.
    dists = [minimum_distance_between_filaments(f_ref, f2) for f2 in f_list]
    if not dists:
        return 1e6
    else:
        return np.min(dists)


# Get random 3D cartesian coordinates
def rand_3D_cartesian_coord(rng, box_size, boundary=None):
    """
    Get a random 3D cartesian point from a uniform distribution.

    Parameters
    ----------
    rng : random number generator
    box_size : array of size (3,)
        the size of the box.

    Returns
    -------
    coord : array of size (3,)
        3D XYZ coordinates.

    """
    coord = np.zeros(box_size.size)
    if boundary is None or boundary['active'] is False:
        for dim in range(coord.size):
            coord[dim] = rng.uniform(0, box_size[dim])
    elif boundary['type'] == 'Tube':
        min_size = boundary['center'] - boundary['radius']
        min_size[2] = 0
        max_size = boundary['center'] + boundary['radius']
        max_size[2] = box_size
        for dim in range(coord.size):
            coord[dim] = rng.uniform(min_size[dim], max_size[dim])
    elif boundary['type'] == 'Sphere':
        min_size = boundary['center'] - boundary['radius']
        max_size = boundary['center'] + boundary['radius']
        for dim in range(coord.size):
            coord[dim] = rng.uniform(min_size[dim], max_size[dim])
    return coord


# Get random 3D orientation vector
def rand_3D_director(rng):
    """
    Get a random 3D unit director from a uniform distribution.

    Parameters
    ----------
    rng : random number generator

    Returns
    -------
    vec : array of size (3,)
        3D unit director.

    """
    # Generate 3 random x,y,z variables. Then normalize.
    vec = rng.normal(size=3)
    vec = vec/np.linalg.norm(vec)
    return vec


# Check if coordinate is within bounds
def check_coordinate_inside_bounds(coord, lower_bound, upper_bound):
    """
    Check if coordinate is within the specified bounds.

    Parameters
    ----------
    coord : array of size (3,)
        3D position to be checked.
    lower_bound : array of size (3,)
    upper_bound : array of size (3,)

    Returns
    -------
    within_bounds : bool
        True is position is within bounds. False if otherwise.

    """
    below_upper = np.all(coord < upper_bound)
    above_lower = np.all(coord > lower_bound)
    within_bounds = below_upper and above_lower
    return within_bounds


# Check if coordinate is inside special boundary
def check_coordinate_inside_boundary(coord, boundary):
    """
    Check if coordinate is within the specified bounds.

    Parameters
    ----------
    coord : array of size (3,)
        3D position to be checked.
    boundary : dict
        fields: active, type, more boundary-specific fields.

    Returns
    -------
    inside : bool
        True is position is within boundary. False if otherwise.

    """
    inside = False
    if boundary['type'] == 'Tube':
        # calc X and Y dist from tube center, ensuring its less than radius
        dist_center = np.linalg.norm(coord[0:2] - boundary['center'][0:2])
        if dist_center < boundary['radius']:
            inside = True

    elif boundary['type'] == 'Sphere':
        # calc X,Y,Z dist from sphere center, ensuring its less than radius
        dist_center = np.linalg.norm(coord - boundary['center'])
        if dist_center < boundary['radius']:
            inside = True

    return inside


# Attempt to create a valid filament
def attempt_create_filament(D, L, f_list, box_size, rng, gid, boundary=None):
    """
    Attempt to add a create a valid filament (inside the sim box).

    Parameters
    ----------
    D : float
        diameter of filament.
    L : float
        length of filament.
    f_list : list
        list of Filament objects to append to.
    gid : int
        Global ID of new filament.
    box_size : array of size (3,)
        size of the sim box.
    rng : random number generator
        random number gen
    boundary: dict
        dict with fields 'active' (bool) and 'type' (string) and some more.

    Returns
    -------
    f_list : list
        List of filament objects.
    gid : int
        New global ID.
    accept_filament : bool
        status of attempt to add filament.

    """
    # Get a 3D random point such that a filament is inside the sim box
    inside = False
    while not inside:

        # random 3D point
        center = rand_3D_cartesian_coord(rng, box_size)

        # random director
        director = rand_3D_director(rng)

        fil = Filament(center, director, L, D/2, gid)

        # Check if inside sim box
        inside = (
            check_coordinate_inside_bounds(
                fil.pos_start, np.array([0, 0, 0]), box_size) and
            check_coordinate_inside_bounds(
                fil.pos_end, np.array([0, 0, 0]), box_size)
            )

        # If boundary active, check if inside boundary
        if (boundary is not None) and boundary['active']:
            inside = (
                inside and
                check_coordinate_inside_boundary(fil.pos_start, boundary) and
                check_coordinate_inside_boundary(fil.pos_end, boundary)
                )

    # find minimum distance to other filaments, and accept if distance is more
    # than the filament diameter
    min_dist = minimum_distance_to_all_filaments(fil, f_list)
    accept_filament = False
    if min_dist > D:
        f_list.append(fil)
        gid += 1
        # print(fil)
        accept_filament = True
    return f_list, gid, accept_filament


# Volume of valid region
def volume_of_valid_region(box_size, boundary=None):
    """
    Calculate the volume of the valid region.

    This can either be the entire simulation volume or the volume inside a
    boundary.

    Parameters
    ----------
    box_size : array of shape (3,)
        dimensions of the sim box.
    boundary : dict, optional
        fields: type, center, radius. The default is None.

    Returns
    -------
    vol: float
        volume of the allowed region.

    """
    if boundary is None or not boundary['active']:
        vol = np.prod(box_size)
    elif boundary['type'] == 'Tube':
        vol = math.pi * boundary['radius']**2 * box_size[-1]
    elif boundary['type'] == 'Sphere':
        vol = (4/3) * math.pi * boundary['radius']**3
    return vol


# =============================================================================
# Code to execute
# =============================================================================

gid = 0  # Begin with a 0 global ID.
f_list = []  # list to store filaments

# If packing fraction is enabled
if use_packing_fraction:

    vol_region = volume_of_valid_region(box_size, boundary)
    vol_enclosed = 0
    current_packing_fraction = 0
    while current_packing_fraction < packing_fraction:
        if sample_len_form_exp_dist:
            L = rand_exponential_length(
                rng, length_mean, length_min, length_max)
        else:
            L = length_mean

        f_list, gid, status = attempt_create_filament(
            diameter_of_filament, L, f_list, box_size, rng, gid, boundary)
        if status:
            vol_enclosed += math.pi * (diameter_of_filament/2)**2 * L
            current_packing_fraction = vol_enclosed/vol_region
            print("", end=f"\rPacking Fraction = {current_packing_fraction}")
            # print('CPF = {}'.format(current_packing_fraction))

else:

    # sample all length values and sort them
    if sample_len_form_exp_dist:
        lens = sorted(
            [rand_exponential_length(rng, length_mean, length_min, length_max)
             for i in range(number_of_filaments)], reverse=True)
    else:
        lens = [length_mean for i in range(number_of_filaments)]

    for idx, el in enumerate(lens):
        status = False  # enables limiting number of fails
        failCount = 0
        while status is False:
            failCount += 1
            f_list, gid, status = attempt_create_filament(
                diameter_of_filament, el, f_list, box_size, rng, gid, boundary)
            if failCount > 5000:
                raise Exception('Failed to add a filament with 5000 trials.\n'
                                'Stopping...\n...\n...\n Total Added Filaments'
                                ' = {}'.format(len(f_list)))

if initialize_proteins:
    p_list = []
    if use_packing_fraction:
        number_of_proteins = len(f_list) * scale_proteins
    for i_protein in range(number_of_proteins):

        outside = True
        while outside:
            coord = rand_3D_cartesian_coord(rng, box_size)
            outside = not check_coordinate_inside_boundary(coord, boundary)
        p_list.append(Protein(coord, gid))
        gid += 1


# Visualization
if do_visualization:
    fig = plt.figure(figsize=(8, 8))
    ax = fig.gca(projection='3d')
    # ax = plt.axes(projection='3d')

    # Plot filaments
    for fil in f_list:
        fil.plot_3D(ax)
    if initialize_proteins:
        for prot in p_list:
            prot.plot_3D(ax)

    ax.set_xlim3d(0, box_size[0])
    ax.set_ylim3d(0, box_size[1])
    ax.set_zlim3d(0, box_size[2])
    if use_packing_fraction:
        figtitle = 'Packing fraction = {}'.format(current_packing_fraction)
    else:
        figtitle = 'Number of filaments = {}'.format(number_of_filaments)
    ax.set(xlabel='X', ylabel='Y', zlabel='Z', title=figtitle)

# Write file
if do_file_write:
    with open(fname, "w") as filer:
        filer.write('# Initial configuration of rods\n#\n')
        for fil in f_list:
            filer.write(fil.get_string_to_write())
    if initialize_proteins:
        with open(pname, "w") as piler:
            piler.write('# Initial configuration of proteins\n#\n')
            for prot in p_list:
                piler.write(prot.get_string_to_write())
