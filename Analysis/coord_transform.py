import numpy as np
from pytransform3d import rotations as pyro
from numba import njit

def transform_coords( coords, origin_new, axis_new):
    # transform coordinates to a new frame of reference with the orivided
    # new origin and the the new basis matrix

    # Translate origin
    coords = coords - origin_new

    # Get vector perp to z-axis and the new z-axis
    vperp = pyro.perpendicular_to_vectors( np.array([0,0,1]), axis_new)

    # Rotate if needed
    if np.linalg.norm( vperp) >0.00001:
        vperp = vperp / np.linalg.norm(vperp)
        theta = pyro.angle_between_vectors( np.array([0,0,1]), axis_new)
        axis_angle = np.concatenate( (vperp, -np.array([theta]) ))
        q = pyro.quaternion_from_axis_angle( axis_angle)

        # transfrom coordinates
        coords = q_prod_coords( coords,q )

    return coords

@njit
def q_prod_coords( coords, q):
    # multiply each row vector in coords by the quaternion q
    for idx in np.arange( coords.shape[0]):
        t = 2 * np.cross(q[1:], coords[idx,:])
        coords[idx,:]+= q[0] * t + np.cross(q[1:], t)
    return coords
