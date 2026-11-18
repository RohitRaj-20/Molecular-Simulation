#!/usr/bin/env python3
"""GROMACS index file for a Martini system built as [nmol CH4C][nion ions][solvent].
Usage: python3 make_ndx_cg.py cg_groups.json index.ndx --nmol 24 --nion 24 --gro system.gro
Groups: System, Solvent, W / DMS (whatever is present), MOL1..MOLn (nmol <= 60), SOLUTE, STER, LINK, ETHER, HEAD,
NPLUS, NITRILE, RING1..RING4 (summed over all copies), ION (bromide), SOLUTE_BR.
In the dimer boxes only the first two molecules exist, so MOL1 and MOL2 are the pull groups."""
import json, argparse
ap = argparse.ArgumentParser()
ap.add_argument('groups'); ap.add_argument('out')
ap.add_argument('--nmol', type=int, required=True); ap.add_argument('--nion', type=int, default=0)
ap.add_argument('--gro', required=True)
a = ap.parse_args()
J = json.load(open(a.groups)); N = J['natoms']; G = J['groups']
def write(f, name, ids):
    f.write(f'[ {name} ]\n')
    for k in range(0, len(ids), 15): f.write(' '.join(str(i) for i in ids[k:k + 15]) + '\n')
gl = open(a.gro).read().splitlines(); n_all = int(gl[1])
res = {}
for k, l in enumerate(gl[2:2 + n_all], 1): res.setdefault(l[5:10].strip(), []).append(k)
solv = [(r, ids) for r, ids in res.items() if r not in (J.get('residue', 'CH4C'), 'ION', 'CH4C')]
with open(a.out, 'w') as f:
    write(f, 'System', list(range(1, n_all + 1)))
    for r, ids in solv: write(f, r, ids)
    write(f, 'Solvent', sorted(i for _, ids in solv for i in ids))
    if a.nmol <= 60:
        for m in range(a.nmol): write(f, f'MOL{m+1}', [m * N + i for i in range(1, N + 1)])
    write(f, 'SOLUTE', list(range(1, a.nmol * N + 1)))
    for name, d in G.items():
        write(f, name, [m * N + i for m in range(a.nmol) for i in d['ids']])
    if a.nion:
        write(f, 'ION', list(range(a.nmol * N + 1, a.nmol * N + a.nion + 1)))
        write(f, 'SOLUTE_BR', list(range(1, a.nmol * N + a.nion + 1)))
print('wrote', a.out)
