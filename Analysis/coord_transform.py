import numpy as np
from pytransform3d import rotations as pyro
from numba import njit

def transform_coords( coords, origin_new, axis_new):
    # transform coordinates to a new frame of reference with the orivided
    # new origin and the the new basis matrix

    world_origin = np.array([0,0,0])
    world_vec = np.identity(3)

    # Translation vector
    T = (world_origin - obj_origin)
    coords = coords + T.transpose()

    # Get vector perp to z-axis and the new z-axis
    vperp = pyro.perpendicular_to_vectors( world_vec[:,-1], axis_new[:,-1])

    # Rotate if needed
    if np.linalg.norm( vperp) >0.00001:
        vperp = vperp / np.linalg.norm(vperp)
        theta = pyro.angle_between_vectors( world_vec[:,-1], obj_vec[:,-1])
        axis_angle = np.concatenate( (vperp, -np.array([theta]) ))
        q = pyro.quaternion_from_axis_angle( axis_angle)

        # check quaternion to ensure its valid
        q = pyro.check_quaternion(q)

        # transfrom coordinates 
        coords = q_prod_coords( coords )

    return coords

@njit
def q_prod_coords( coords, q):
    # multiply each row vector in coords by the quaternion q
    for idx in np.arange( coords.shape[0]):
        t = 2 * np.cross(q[1:], coords[idx,:])
        coords[idx,:]+= q[0] * t + np.cross(q[1:], t)
    return coords
