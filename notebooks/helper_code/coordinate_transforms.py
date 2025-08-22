import numpy as np

def xyz_to_rz(vec):
    """ Convert a 3-vector to cylindrical coordinates (r, z)."""
    r = r = np.linalg.norm(vec[:2]) 
    z = vec[2]
    return np.array([r, z])


def rz_to_xyz(vec, normal):
    """Convert a 2-vector to cartesian coordinates. (x, y, z)
    """
    r = vec[0]
    z = vec[1]

    xy_unit = normal[:2] / np.linalg.norm(normal[:2]) if np.linalg.norm(normal[:2]) > 0 else np.array([1, 0])
    return np.array([r * xy_unit[0], r * xy_unit[1], z])
