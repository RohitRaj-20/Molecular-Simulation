#!/usr/bin/env python3
"""Two Martini molecules in a cubic box for the dimer PMF (CG version of make_dimer.py).
Usage:
  python3 make_dimer_cg.py CH4C_cg.gro CH4C.itp solutes.gro --mode parallel --sep 4.5 --box 10
Modes (the separation axis is rotated onto z, the pulling direction):
  parallel      copy 2 is a translated copy of copy 1 along the short (face-to-face) axis
  antiparallel  copy 2 is rotated 180 degrees about that axis (long axes reversed)
  endtoend      copy 2 is a translated copy along the long axis
Prints the COM distance, the closest bead contact and the two pull-group pbcatom numbers."""
import argparse, numpy as np
from cg_common import parse_itp_atoms, read_gro

def build(xyz, mass, mode, sep):
    w = mass if mass.sum() > 0 else np.ones(len(xyz))
    com = (xyz * w[:, None]).sum(0) / w.sum(); X = xyz - com
    ev, v = np.linalg.eigh(np.cov(X.T))
    e3, e2, e1 = v[:, 0], v[:, 1], v[:, 2]
    zax, xax = (e1, e2) if mode == 'endtoend' else (e3, e1)
    yax = np.cross(zax, xax); R = np.vstack([xax, yax, zax])
    assert abs(np.linalg.det(R) - 1) < 1e-9
    A = X @ R.T
    B = A * np.array([-1, -1, 1]) if mode == 'antiparallel' else A.copy()
    return A + np.array([0, 0, -sep / 2]), B + np.array([0, 0, sep / 2]), w

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('gro'); ap.add_argument('itp'); ap.add_argument('out')
    ap.add_argument('--mode', choices=['parallel', 'antiparallel', 'endtoend'], default='parallel')
    ap.add_argument('--sep', type=float, default=4.5); ap.add_argument('--box', type=float, default=10.0)
    a = ap.parse_args()
    names, resn, xyz, _ = read_gro(a.gro); _, _, _, mass = parse_itp_atoms(a.itp)
    A, B, w = build(xyz, mass, a.mode, a.sep)
    n = len(names); res = resn[0]
    with open(a.out, 'w') as f:
        f.write(f'dimer\n{2 * n}\n'); k = 0
        for r, P in enumerate((A, B), 1):
            for i in range(n):
                k += 1; x, y, z = P[i] + a.box / 2
                f.write(f'{r:5d}{res:<5s}{names[i]:>5s}{k % 100000:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n')
        f.write(f'{a.box:10.5f}{a.box:10.5f}{a.box:10.5f}\n')
    comA = (A * w[:, None]).sum(0) / w.sum(); comB = (B * w[:, None]).sum(0) / w.sum()
    d = np.linalg.norm(A[:, None] - B[None], axis=2)
    cA = int(np.argmin(np.where(w > 0, np.linalg.norm(A - comA, axis=1), 1e9))) + 1
    print(f'mode {a.mode}: COM distance {np.linalg.norm(comB - comA):.3f} nm, closest bead contact {d.min():.3f} nm, '
          f'extent {(A.max(0) - A.min(0)).round(2)} nm')
    print(f'pull-group1-pbcatom = {cA}\npull-group2-pbcatom = {n + cA}')
    if a.sep > a.box / 2 * 0.98: print('WARNING: sep is larger than half the box')
