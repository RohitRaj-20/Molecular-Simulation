"""Geometry helpers for the aggregation analysis (numpy + scipy only, cubic or orthorhombic boxes)."""
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

def molecule_contacts(heavy_xyz, mol_of_atom, box, cutoff):
    """Set of (i, j) molecule pairs with any heavy-atom pair closer than cutoff (periodic)."""
    tree = cKDTree(np.mod(heavy_xyz, box), boxsize=box)
    pairs = tree.query_pairs(cutoff, output_type='ndarray')
    if len(pairs) == 0: return np.zeros((0, 2), int)
    mi, mj = mol_of_atom[pairs[:, 0]], mol_of_atom[pairs[:, 1]]
    keep = mi != mj
    mp = np.sort(np.c_[mi[keep], mj[keep]], axis=1)
    return np.unique(mp, axis=0)

def clusters_from_pairs(nmol, mol_pairs):
    if len(mol_pairs) == 0: return np.arange(nmol)
    A = coo_matrix((np.ones(len(mol_pairs)), (mol_pairs[:, 0], mol_pairs[:, 1])), shape=(nmol, nmol))
    return connected_components(A, directed=False)[1]

def unwrap_cluster(mol_xyz, members, mol_pairs, box):
    """mol_xyz: list of (natoms,3) arrays of WHOLE molecules. Returns the members shifted by
    box vectors so the cluster is contiguous. Valid while neighbour COM distance < box/2."""
    members = list(members); mset = set(members)
    nbr = {m: [] for m in members}
    for i, j in mol_pairs:
        if i in mset and j in mset: nbr[i].append(j); nbr[j].append(i)
    out = {members[0]: mol_xyz[members[0]].copy()}
    stack = [members[0]]
    while stack:
        u = stack.pop(); cu = out[u].mean(0)
        for v in nbr[u]:
            if v in out: continue
            shift = -box * np.round((mol_xyz[v].mean(0) - cu) / box)
            out[v] = mol_xyz[v] + shift; stack.append(v)
    return out

def gyration(xyz, mass):
    com = (xyz * mass[:, None]).sum(0) / mass.sum()
    d = xyz - com
    T = (d[:, :, None] * d[:, None, :] * mass[:, None, None]).sum(0) / mass.sum()
    lam = np.sort(np.linalg.eigvalsh(T))[::-1]
    rg = np.sqrt(lam.sum())
    kappa2 = 1.5 * (lam ** 2).sum() / lam.sum() ** 2 - 0.5 if lam.sum() > 0 else 0.0
    return rg, kappa2

def ring_geometry(ring_xyz):
    c = ring_xyz.mean(0)
    n = np.linalg.svd(ring_xyz - c)[2][-1]
    return c, n

def stacked_pairs(centroids, normals, mol_ids, box, d_cut=0.55, ang_cut=30.0):
    """Number of intermolecular ring pairs closer than d_cut with normals within ang_cut of parallel."""
    tree = cKDTree(np.mod(centroids, box), boxsize=box)
    cnt, partners = 0, set()
    for i, j in tree.query_pairs(d_cut):
        if mol_ids[i] == mol_ids[j]: continue
        cosang = abs(float(np.dot(normals[i], normals[j])))
        if np.degrees(np.arccos(min(1.0, cosang))) <= ang_cut:
            cnt += 1; partners.add(mol_ids[i]); partners.add(mol_ids[j])
    return cnt, len(partners)

def ion_pairs(cation_xyz, ion_xyz, box, cut=0.5):
    if len(ion_xyz) == 0: return 0.0
    tree = cKDTree(np.mod(ion_xyz, box), boxsize=box)
    return float(np.mean([len(x) for x in tree.query_ball_point(np.mod(cation_xyz, box), cut)]))
