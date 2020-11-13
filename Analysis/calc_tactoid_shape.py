import numpy as np
import open3d as o3d

def calc_aspect_ratio(coords):
    # Calculate aspect ratio of bounding box of a collection of 3D points

    if not coords.any():
        return np.nan

    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector( coords)
    obb = pcd.get_oriented_bounding_box()
    xx = obb.extent
    leng = max(xx)
    width = np.mean(sorted(xx)[0:2])
    aspect_ratio = leng / width
    return aspect_ratio, leng, width
