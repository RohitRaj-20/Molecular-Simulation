#!/usr/bin/env python3
"""Compare the CG bonded terms with an all-atom trajectory of the same molecule and propose new values.
Usage:
  gmx trjconv -s aa_md.tpr -f aa_md.xtc -n index.ndx -pbc mol -o aa_solute.gro     # group: the molecule only (130 atoms)
  python3 cg_bonded_from_aa.py --mol2 Chol4CSN_clean.mol2 --map chol4csn_cg_map.json --index cg_index.json \
        --gro aa_solute.gro --out bonded_overrides.json
Every frame is mapped to the head/linker beads (centre of the mapped heavy atoms, rings NOT scaled), and the
distributions of every bond, angle and ring-ring torsion are reported. Proposed values:
  bonds    keep the CG length (the Martini ring size makes AA-mapped lengths too short); k = kT / sigma^2, clipped to 3000..20000
  angles   theta0 = AA circular mean; k = kT / (sin(theta0) sigma)^2, clipped to 20..400   (angle function 2)
  dihedral for each ring-ring bond, if the AA distribution is unimodal (circular sd < 30 deg): function 2, k = kT/sigma^2 clipped 5..200
Feed the result back:  python3 cg_build.py mol2 map.json --outdir . --overrides bonded_overrides.json
Only numpy is needed."""
import json, argparse
import numpy as np
from cg_common import gro_frames
from cg_build import read_mol2, dihedral, angle, KT

ap = argparse.ArgumentParser()
ap.add_argument('--mol2', required=True); ap.add_argument('--map', required=True); ap.add_argument('--index', required=True)
ap.add_argument('--gro', required=True); ap.add_argument('--out', default='bonded_overrides.json')
ap.add_argument('--stride', type=int, default=1)
a = ap.parse_args()
names, _, _ = read_mol2(a.mol2); idx = {n: i for i, n in enumerate(names)}; NA = len(names)
I = json.load(open(a.index)); order = I['order']; nid = {n: i + 1 for i, n in enumerate(order)}
M = json.load(open(a.map))
beads = [b for b in I['beads']]
ids = {b['name']: [idx[n] for n in b['atoms']] for b in beads}
ids['ETH'] = [idx[n] for n in M['steroid']['anchor_head']]      # ETH bead sits on C3 + ether O
ring_of = {b['name']: b['ring'] for b in beads if 'ring' in b}
def mate(b, other):          # ring neighbours of bead b that are not 'other'
    return [x for x in ring_of if ring_of[x] == ring_of[b] and x != b and x != other]
B = [(order[b['i'] - 1], order[b['j'] - 1]) for b in I['bonds']]
A = [(order[x['i'] - 1], order[x['j'] - 1], order[x['k'] - 1]) for x in I['angles']]
D = []
_x0 = next(gro_frames(a.gro))[1][:NA]
_P0 = {b: _x0[ids[b]].mean(0) for b in ids}
def best_mate(p, q):         # ring neighbour of p that is NOT nearly collinear with p-q (a collinear triple has no torsion)
    return max(mate(p, q), key=lambda m: abs(np.sin(np.radians(angle(_P0[m], _P0[p], _P0[q])))))
for p, q in B:
    if p in ring_of and q in ring_of and ring_of[p] != ring_of[q]:
        D.append((best_mate(p, q), p, q, best_mate(q, p)))
S = {('b',) + x: [] for x in B}; S.update({('a',) + x: [] for x in A}); S.update({('d',) + x: [] for x in D})
nf = 0
for k, (t, xyz, box) in enumerate(gro_frames(a.gro)):
    if k % a.stride: continue
    nmol = len(xyz) // NA
    for m in range(nmol):
        X = xyz[m * NA:(m + 1) * NA]
        P = {b: X[ids[b]].mean(0) for b in ids}
        for x in B: S[('b',) + x].append(np.linalg.norm(P[x[0]] - P[x[1]]))
        for x in A: S[('a',) + x].append(angle(P[x[0]], P[x[1]], P[x[2]]))
        for x in D: S[('d',) + x].append(dihedral(*[P[y] for y in x]))
    nf += 1
print(f'{nf} frames')
ov = dict(bonds=[], angles=[], dihedrals=[])
print('\nBONDS (AA-mapped mean +/- sd, nm)   current CG length is kept')
for x in B:
    v = np.array(S[('b',) + x]); sd = v.std()
    kk = float(np.clip(KT / max(sd, 1e-3) ** 2, 3000, 20000))
    print(f'  {x[0]:5s}-{x[1]:5s} {v.mean():.3f} +/- {sd:.3f}   proposed k = {kk:.0f}')
    cur = next(b for b in I['bonds'] if (order[b['i'] - 1], order[b['j'] - 1]) == x)
    ov['bonds'].append(dict(i=cur['i'], j=cur['j'], length=cur['length'], k=kk))
print('\nANGLES (deg)')
for x in A:
    v = np.radians(np.array(S[('a',) + x])); mean = np.degrees(np.arctan2(np.sin(v).mean(), np.cos(v).mean())); sd = v.std()
    kk = float(np.clip(KT / (np.sin(np.radians(mean)) * max(sd, 1e-2)) ** 2, 20, 400))
    print(f'  {x[0]:5s}{x[1]:5s}{x[2]:5s} {mean:7.1f} +/- {np.degrees(sd):5.1f}   proposed k = {kk:.0f}')
    ov['angles'].append(dict(i=nid[x[0]], j=nid[x[1]], k=nid[x[2]], theta=float(mean), k_force=kk))
print('\nRING-RING TORSIONS (deg)')
for x in D:
    v = np.radians(np.array(S[('d',) + x])); R = np.hypot(np.sin(v).mean(), np.cos(v).mean())
    mean = np.degrees(np.arctan2(np.sin(v).mean(), np.cos(v).mean())); csd = np.degrees(np.sqrt(-2 * np.log(max(R, 1e-9))))
    msg = f'  {x[0]}-{x[1]}-{x[2]}-{x[3]}  mean {mean:7.1f}  circular sd {csd:5.1f}'
    if csd < 30:
        kk = float(np.clip(KT / np.radians(csd) ** 2, 5, 200)); msg += f'   proposed k = {kk:.0f}'
        ov['dihedrals'].append(dict(i=nid[x[0]], j=nid[x[1]], k=nid[x[2]], l=nid[x[3]], phi=float(mean), k_force=kk,
                                    note='ring-ring torsion from AA'))
    else:
        msg += '   multimodal or flat: no torsion term written, inspect the histogram'
    print(msg)
json.dump(ov, open(a.out, 'w'), indent=1); print('\nwrote', a.out)
