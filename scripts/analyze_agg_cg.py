#!/usr/bin/env python3
"""Per-frame aggregation metrics for the Martini N-molecule box (numpy + scipy; MDAnalysis only for .xtc input).
Prepare the solute-only trajectory (molecules whole):
  gmx trjconv -s md.tpr -f md.xtc -n index.ndx -pbc mol -o solute_whole.xtc       # group SOLUTE_BR
  gmx trjconv -s md.tpr -f md.xtc -n index.ndx -pbc mol -dump 0 -o solute.gro      # same group
Run:
  python3 analyze_agg_cg.py --top solute.gro --traj solute_whole.xtc --groups cg_groups.json --itp CH4C.itp \
      --nmol 24 --nion 24 --cutoff 0.6 --out agg_rep1.csv
A .gro multi-frame file (gmx trjconv -o solute_whole.gro) is read natively, with no MDAnalysis needed.
Atom order must be [nmol x CH4C beads][nion bromides].
Cluster rule: two molecules share a cluster when any two of their beads are closer than --cutoff (nm).
Ring stacking: two ring triplets of different molecules whose centroids are closer than --ring-cut and whose
normals are within 30 degrees of parallel.  Ion pairing: mean number of Br- within --ion-cut of each N+ bead.
Columns: time_ps, n_clusters, n_monomers, largest, mean_size_w, rg_largest_nm, kappa2_largest, stacked_pairs,
mol_in_stacks, br_per_cation  (same as the all-atom analysis, so summarize_reps.py works unchanged)."""
import argparse, json, csv
import numpy as np
from agg_core import *
from cg_common import parse_itp_atoms, gro_frames

def frames(a):
    if a.traj.endswith('.gro'):
        for k, (t, xyz, box) in enumerate(gro_frames(a.traj)):
            if k % a.stride == 0: yield t, xyz, box
    else:
        import MDAnalysis as mda
        u = mda.Universe(a.top, a.traj)
        for ts in u.trajectory[::a.stride]: yield ts.time, u.atoms.positions / 10.0, ts.dimensions[:3] / 10.0

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--top'); ap.add_argument('--traj', required=True)
    ap.add_argument('--groups', required=True); ap.add_argument('--itp', required=True)
    ap.add_argument('--nmol', type=int, required=True); ap.add_argument('--nion', type=int, default=0)
    ap.add_argument('--cutoff', type=float, default=0.6)
    ap.add_argument('--ring-cut', type=float, default=0.6); ap.add_argument('--ion-cut', type=float, default=0.6)
    ap.add_argument('--stride', type=int, default=1); ap.add_argument('--out', default='agg.csv')
    a = ap.parse_args()
    G = json.load(open(a.groups)); N = G['natoms']; g = G['groups']
    mass1 = parse_itp_atoms(a.itp)[3]
    assert len(mass1) == N, 'itp and groups.json disagree on the number of beads'
    mass1 = np.where(mass1 > 0, mass1, 1.0)               # virtual sites carry no mass in the topology; keep them in Rg
    ring_ids = [np.array(g[k]['ids']) - 1 for k in g if k.startswith('RING')]
    nplus = g['NPLUS']['ids'][0] - 1
    mol_of = np.repeat(np.arange(a.nmol), N)
    out = open(a.out, 'w', newline=''); w = csv.writer(out)
    w.writerow(['time_ps', 'n_clusters', 'n_monomers', 'largest', 'mean_size_w', 'rg_largest_nm',
                'kappa2_largest', 'stacked_pairs', 'mol_in_stacks', 'br_per_cation'])
    for t, pos, box in frames(a):
        assert len(pos) >= a.nmol * N + a.nion, 'trajectory has fewer beads than nmol*N+nion'
        sol = pos[:a.nmol * N].reshape(a.nmol, N, 3)
        pairs = molecule_contacts(sol.reshape(-1, 3), mol_of, box, a.cutoff)
        lab = clusters_from_pairs(a.nmol, pairs); sizes = np.bincount(lab)
        big = int(np.argmax(sizes)); members = np.where(lab == big)[0]
        un = unwrap_cluster([sol[m] for m in range(a.nmol)], members, pairs, box)
        X = np.vstack([un[m] for m in sorted(un)]); Mm = np.tile(mass1, len(un))
        rg, k2 = gyration(X, Mm) if len(un) > 1 else (np.nan, np.nan)
        cen, nor, ids = [], [], []
        for m in range(a.nmol):
            for r in ring_ids:
                c, nv = ring_geometry(sol[m][r]); cen.append(c); nor.append(nv); ids.append(m)
        sp, sm = stacked_pairs(np.array(cen), np.array(nor), np.array(ids), box, d_cut=a.ring_cut)
        ions = pos[a.nmol * N:a.nmol * N + a.nion]
        w.writerow([f'{t:.1f}', len(sizes), int((sizes == 1).sum()), int(sizes.max()),
                    f'{(sizes ** 2).sum() / a.nmol:.3f}', f'{rg:.3f}', f'{k2:.3f}', sp, sm,
                    f'{ion_pairs(sol[:, nplus, :], ions, box, cut=a.ion_cut):.3f}'])
    out.close(); print('wrote', a.out)

if __name__ == '__main__':
    main()
