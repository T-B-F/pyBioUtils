import numpy as np

def compute_rmsd(coords1, coords2, return_struct=False):
    if not isinstance(coords1, np.ndarray) or not isinstance(coords2, np.ndarray):
        raise ValueError("coordinates must be numpy array")
    if not coords1.shape == coords2.shape:
        raise ValueError("coordinates must be have same shape")
    if len(coords1.shape) != 2:
        raise ValueError("coordinates must be two dimensionnal")

    # Kabsch algo
    cp_coords1 = coords1.copy()
    cp_coords2 = coords2.copy()
    # center 
    cp_coords1 -= cp_coords1.sum(axis=0) / cp_coords1.shape[0]
    cp_coords2 -= cp_coords2.sum(axis=0) / cp_coords2.shape[0]
    V, S, Wt = np.linalg.svd(np.dot(cp_coords1.T, cp_coords2))
    reflect = np.linalg.det(V) * np.linalg.det(Wt)
    if reflect < -1.0:
        S[-1] = -S[-1]
        V[:,-1] = -V[:,-1]

    E0 = ((cp_coords2)**2).sum() + ((cp_coords1)**2).sum()
    rmsd = E0 - (2.0 * sum(S))
    rmsd = np.sqrt(abs(rmsd / cp_coords2.shape[0]))
    # rotate, translate, center
    U = np.dot(V, Wt)
    cp_coords1 = np.dot(cp_coords1 , U)
    if return_struct:
        return rmsd, cp_coords1, cp_coords2
    else:
        return rmsd

def compute_rmsf(coords, target):
    """ compute fluctuation, trajectories should already been superposed
    """
    idx = np.arange(target.shape[0])
    return np.sqrt(3*np.mean((coords[:, idx, :] - target) ** 2, axis=(0, 2)))
    
