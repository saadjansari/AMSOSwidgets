import numpy as np
import open3d as o3d

def calc_aspect_ratio( coords):
    # Calculate aspect ratio of bounding box of a collection of 3D points

    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector( coords)
    obb = pcd.get_oriented_bounding_box()

    xx = obb.extent
    aspect_ratio = max( xx) / np.mean(sorted(xx)[0:2])
    return aspect_ratio

# def remove_outliers_point_cloud( point_cloud):
    # remove outliers from the point cloud

# @njit
# def get_linspace_start_end( c0, c1):
    # calculate a collection of points lying between c0 and c1

