import numpy as np
import open3d as o3d


def calc_aspect_ratio(coords):
    # Calculate aspect ratio of bounding box of a collection of 3D points

    if not coords:
        return np.nan
    xx = obb.extent
    aspect_ratio = max( xx) / np.mean(sorted(xx)[0:2])
    return aspect_ratio
