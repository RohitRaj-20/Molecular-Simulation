#!/usr/bin/env python3
"""Lattice-built Martini solvent box: water (W beads, 4 waters each), DMSO (2 beads, 1 molecule) or a mixture.
Usage:
  python3 make_solvent_cg.py --n 10 --x 0.0  --out water_raw.gro     # pure water
  python3 make_solvent_cg.py --n 10 --x 1.0  --out dmso_raw.gro      # pure DMSO
  python3 make_solvent_cg.py --n 10 --x 0.3  --out mix30_raw.gro     # X(DMSO) = 0.3 of all (real) water+DMSO molecules
The box has n^3 sites (default spacing 0.49 nm). x is the mole fraction counting one W bead as 4 water molecules:
  n_DMSO = 4 x N / (1 + 3 x),  n_W = N - n_DMSO.
A W bead (4 waters, ~0.120 nm^3) and a DMSO molecule (~0.118 nm^3) occupy nearly the same volume, so the NPT run
only has to adjust the box slightly. Prints the counts and the topology [ molecules ] block."""
import argparse, numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('--n', type=int, default=10); ap.add_argument('--x', type=float, required=True)
ap.add_argument('--spacing', type=float, default=0.49); ap.add_argument('--seed', type=int, default=1)
ap.add_argument('--out', required=True)
a = ap.parse_args()
rng = np.random.default_rng(a.seed)
N = a.n ** 3; nD = int(round(4 * a.x * N / (1 + 3 * a.x))); nW = N - nD
sites = np.array([[i, j, k] for i in range(a.n) for j in range(a.n) for k in range(a.n)], float) * a.spacing + a.spacing / 2
rng.shuffle(sites)
box = a.n * a.spacing
lines = []; atom = 0; res = 0
for s in sites[:nW]:
    res += 1; atom += 1
    lines.append(f'{res % 100000:5d}{"W":<5s}{"W":>5s}{atom % 100000:5d}{s[0]:8.3f}{s[1]:8.3f}{s[2]:8.3f}')
for s in sites[nW:]:
    res += 1
    u = rng.normal(size=3); u /= np.linalg.norm(u)
    for nm, p in (('S1', s), ('O2', s + 0.30 * u)):
        atom += 1; p = np.mod(p, box)
        lines.append(f'{res % 100000:5d}{"DMS":<5s}{nm:>5s}{atom % 100000:5d}{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}')
with open(a.out, 'w') as f:
    f.write(f'Martini solvent x_DMSO={a.x}\n{atom}\n' + '\n'.join(lines) + f'\n{box:10.5f}{box:10.5f}{box:10.5f}\n')
nwat = 4 * nW
print(f'box {box:.3f} nm: {nW} W beads (= {nwat} waters), {nD} DMSO; X(DMSO) = {nD / (nD + nwat):.3f}')
print('[ molecules ]'); print(f'W    {nW}') if nW else None; print(f'DMSO {nD}') if nD else None
