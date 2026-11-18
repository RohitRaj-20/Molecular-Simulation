#!/usr/bin/env python3
"""Build the Martini 3 model of Chol4CSN (residue CH4C) from the all-atom mol2.

Usage:
  python3 cg_build.py ../params/Chol4CSN_clean.mol2 ../params/chol4csn_cg_map.json --outdir .
  python3 cg_build.py mol2 map.json --outdir . --overrides bonded_overrides.json   # after cg_bonded_from_aa.py

Writes in --outdir:
  CH4C.itp          Martini 3 topology (cholesterol part from the published Martini 3 CHOL model, ROH bead
                    retyped as an ether bead and bonded to the linker; head group mapped from your mol2)
  CH4C_cg.gro       one molecule, coordinates in nm, centred at the origin
  cg_groups.json    bead groups for make_ndx_cg.py / analyze_agg_cg.py
  cg_index.json     bead, bond and angle lists used by cg_bonded_from_aa.py
Only numpy is needed.

The steroid bead positions are CHOL template positions fitted onto your all-atom steroid.
The CHOL parameters are copied from martini_v3.0.0_sterols_v1.itp (Borges-Araujo et al. 2023, Apache 2.0).
"""
import json, argparse, os
import numpy as np

# ---- Martini 3 cholesterol (CHOL) as published; beads 1..9 -----------------------------------------
CHOL_NAMES = ['ETH', 'R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'C1', 'C2']      # bead 1 is ROH in the original
CHOL_TYPES = ['N2a', 'SC4', 'SC3', 'SC3', 'SC3', 'TC2', 'TC2', 'C2', 'C2']  # bead 1 type overwritten by the map
CHOL_MASS  = [0.0, 159.8, 157.0, 0.0, 0.0, 0.0, 0.0, 115.2, 72.0]
CHOL_CONSTR = [(8, 3, 0.75012), (8, 2, 0.78504), (3, 2, 0.34797)]         # (i, j, length nm), 1-based
CHOL_BOND = (8, 9, 0.440, 15000)
CHOL_ANGLE = (5, 4, 9, 1, 99.0, 250)
CHOL_DIH = (7, 5, 4, 9, 2, -70.0, 50)
CHOL_VS = [(1, 1.08999, 0.35891, 0.22947), (4, -0.41067, 0.83597, -0.05311), (5, 0.69138, -0.10352, 0.17476),
           (6, 0.77524, 0.22413, 0.92700), (7, 0.06219, 0.32398, 0.76967)]  # constructed from beads 8, 3, 2 (3out)
KT = 0.0083144626 * 300.0

def read_mol2(path):
    t = open(path).read().replace('\r', '')
    rows = [l.split() for l in t.split('@<TRIPOS>ATOM')[1].split('@<TRIPOS>BOND')[0].strip().splitlines()]
    names = [r[1] for r in rows]
    xyz = np.array([[float(v) for v in r[2:5]] for r in rows]) / 10.0
    el = [r[5].split('.')[0] for r in rows]
    return names, xyz, el

def chol_template():
    """Coordinates (nm) of the nine CHOL beads from the three real beads 8, 3, 2 and the 3out virtual sites."""
    P = np.zeros((9, 3))
    a, b, c = 0.75012, 0.78504, 0.34797
    x = (b * b + a * a - c * c) / (2 * a); y = np.sqrt(b * b - x * x)
    P[7] = [0, 0, 0]; P[2] = [a, 0, 0]; P[1] = [x, y, 0]
    rij, rik = P[2] - P[7], P[1] - P[7]
    for site, A, B, C in CHOL_VS:
        P[site - 1] = P[7] + A * rij + B * rik + C * np.cross(rij, rik)
    ring = P[1:7].mean(0)
    d = P[7] - ring; d /= np.linalg.norm(d)
    P[8] = P[7] + CHOL_BOND[2] * d
    return P

def rot_about(axis, ang):
    axis = axis / np.linalg.norm(axis); K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(ang) * K + (1 - np.cos(ang)) * K @ K

def rot_a_to_b(a, b):
    a = a / np.linalg.norm(a); b = b / np.linalg.norm(b)
    v = np.cross(a, b); s = np.linalg.norm(v); c = float(a @ b)
    if s < 1e-9: return np.eye(3) if c > 0 else rot_about(np.cross(a, [1, 0, 0]) if abs(a[0]) < .9 else np.cross(a, [0, 1, 0]), np.pi)
    K = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + K + K @ K * ((1 - c) / s ** 2)

def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - (b0 @ b1n) * b1n; w = b2 - (b2 @ b1n) * b1n
    x = v @ w; y = np.cross(b1n, v) @ w
    return float(np.degrees(np.arctan2(y, x)))

def angle(p, q, r):
    u, v = p - q, r - q
    return np.degrees(np.arccos(np.clip(u @ v / np.linalg.norm(u) / np.linalg.norm(v), -1, 1)))

def relax(X, bonds, constr, angles, nbr, n_fixed=9, rmin=0.34, steps=4000):
    """Pull the head/linker beads to their bond, constraint and angle targets and push non-bonded beads
    (graph distance >= 2) apart. Steroid beads 1..n_fixed stay where they are. Gives a sane start geometry only."""
    from scipy.optimize import minimize
    n = len(X); free = np.arange(n_fixed, n)
    X0 = X.copy()
    # graph distance <= 1 pairs are excluded from the repulsion
    near = {(b['i'] - 1, b['j'] - 1) for b in bonds} | {(c['i'] - 1, c['j'] - 1) for c in constr}
    near |= {(j, i) for (i, j) in near}
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n) if (i, j) not in near and not (i < n_fixed and j < n_fixed)]
    P = np.array(pairs)
    def energy(v):
        Y = X0.copy(); Y[free] = v.reshape(-1, 3)
        E = 0.0
        for b in bonds:
            E += 0.5 * 1000 * (np.linalg.norm(Y[b['i'] - 1] - Y[b['j'] - 1]) - b['length']) ** 2
        for c in constr:
            E += 0.5 * 5000 * (np.linalg.norm(Y[c['i'] - 1] - Y[c['j'] - 1]) - c['length']) ** 2
        for a in angles:
            u, w = Y[a['i'] - 1] - Y[a['j'] - 1], Y[a['k'] - 1] - Y[a['j'] - 1]
            ct = u @ w / np.linalg.norm(u) / np.linalg.norm(w)
            E += 0.5 * 20 * (ct - np.cos(np.radians(a['theta']))) ** 2 * 10
        d = np.linalg.norm(Y[P[:, 0]] - Y[P[:, 1]], axis=1)
        E += 0.5 * 2000 * np.sum(np.clip(rmin - d, 0, None) ** 2)
        return E
    r = minimize(energy, X0[free].ravel(), method='L-BFGS-B', options=dict(maxiter=steps, maxfun=10 * steps))
    X0[free] = r.x.reshape(-1, 3)
    return X0

def build(mol2, M, overrides=None):
    names, xyz, el = read_mol2(mol2)
    idx = {n: i for i, n in enumerate(names)}
    A = lambda lst: xyz[[idx[n] for n in lst]]
    cog = lambda lst: A(lst).mean(0)
    st = M['steroid']
    # ---- head / linker beads from the all-atom geometry, rings scaled to the Martini ring size ----
    beads = list(M['beads']); pos = {}
    rings = {}
    for b in beads:
        pos[b['name']] = cog(b['atoms'])
        if 'ring' in b: rings.setdefault(b['ring'], []).append(b['name'])
    pos0 = {k: v.copy() for k, v in pos.items()}          # unscaled all-atom mapped positions (angle targets)
    scale = {}
    for r, bl in rings.items():
        assert len(bl) == 3, r
        ring_atoms = sorted(set(a for b in beads if b.get('ring') == r for a in b['atoms']))
        cen = A(ring_atoms).mean(0)
        side = np.mean([np.linalg.norm(pos[bl[i]] - pos[bl[(i + 1) % 3]]) for i in range(3)])
        scale[r] = M['ring_bond_nm'] / side
        for n in bl: pos[n] = cen + scale[r] * (pos[n] - cen)
    # ---- steroid: template fitted on the all-atom steroid ----
    T = chol_template()
    q1 = cog(st['anchor_head']); q2 = cog(st['anchor_tail'])
    p1 = T[0]; p2 = T[8]
    R0 = rot_a_to_b(p2 - p1, q2 - q1)
    ster_xyz = A(st['heavy_atoms'])
    best = None
    for roll in np.radians(np.arange(0, 360, 5)):
        R = rot_about(q2 - q1, roll) @ R0
        X = (T - p1) @ R.T + q1
        d = np.linalg.norm(X[1:8, None] - ster_xyz[None], axis=2).min(1).sum()
        if best is None or d < best[0]: best = (d, X)
    S = best[1]                                           # 9 steroid beads
    pos['ETH'] = S[0]; pos0['ETH'] = S[0]
    # ---- assemble ----
    order = list(CHOL_NAMES) + [b['name'] for b in beads]
    types = list(CHOL_TYPES); types[0] = st['eth_type']
    types += [b['type'] for b in beads]
    charges = [0.0] * 9 + [float(b['charge']) for b in beads]
    mass = list(CHOL_MASS) + [None] * len(beads)
    X = np.vstack([S] + [pos[b['name']][None] for b in beads])
    X = X - X.mean(0)
    nid = {n: i + 1 for i, n in enumerate(order)}
    # ---- bonds, constraints, angles ----
    bonds = []
    lmin = M.get('link_min_nm', 0.0)
    for a_, b_, k in M['bonds']:
        i, j = nid[a_], nid[b_]
        bonds.append(dict(i=i, j=j, length=float(max(lmin, np.linalg.norm(X[i - 1] - X[j - 1]))), k=float(k)))
    constr = []
    for r, bl in rings.items():
        for u, v in ((0, 1), (1, 2), (0, 2)):
            constr.append(dict(i=nid[bl[u]], j=nid[bl[v]], length=M['ring_bond_nm']))
    nbr = {i: [] for i in range(1, len(order) + 1)}
    for b in bonds: nbr[b['i']].append(b['j']); nbr[b['j']].append(b['i'])
    for c in constr: nbr[c['i']].append(c['j']); nbr[c['j']].append(c['i'])
    ring_of = {nid[b['name']]: b['ring'] for b in beads if 'ring' in b}
    kdef = M['default_angle_k']
    angles = []
    for c in range(10, len(order) + 1):                    # beads 10.. (steroid has its own CHOL terms; ETH has no angle)
        ns = sorted(nbr[c])
        for u in range(len(ns)):
            for v in range(u + 1, len(ns)):
                i, k = ns[u], ns[v]
                if c in ring_of and ring_of.get(i) == ring_of[c] and ring_of.get(k) == ring_of[c]: continue
                kind = 'branch' if len(nbr[c]) >= 3 and c not in ring_of else ('ring_exit' if c in ring_of else 'chain')
                P0 = lambda t: pos0[order[t - 1]] if order[t - 1] in pos0 else X[t - 1]
                angles.append(dict(i=i, j=c, k=k, theta=float(angle(P0(i), P0(c), P0(k))), k_force=float(kdef[kind])))
    if overrides:
        for o in overrides.get('bonds', []):
            for b in bonds:
                if {b['i'], b['j']} == {o['i'], o['j']}: b.update(length=o['length'], k=o['k'])
        for o in overrides.get('angles', []):
            for a in angles:
                if a['j'] == o['j'] and {a['i'], a['k']} == {o['i'], o['k']}: a.update(theta=o['theta'], k_force=o['k'])
    dihedrals = list(overrides.get('dihedrals', [])) if overrides else []
    X = relax(X, bonds, constr, angles, nbr, n_fixed=9)
    # keep the E/Z configuration of the C=C bead: one improper-type dihedral over its three neighbours,
    # reference value taken from the all-atom conformer (unscaled mapped positions)
    if all(n in nid for n in ('P1_C', 'CN', 'VIN', 'P2_A')):
        q = [nid[n] for n in ('P1_C', 'CN', 'VIN', 'P2_A')]
        phi0 = dihedral(*[pos0[order[t - 1]] for t in q])
        dihedrals.insert(0, dict(i=q[0], j=q[1], k=q[2], l=q[3], phi=phi0, k_force=50.0,
                                 note='keeps the E/Z configuration at the vinyl bead'))
    groups = {'STER': list(range(1, 10)), 'LINK': [nid['LK']], 'ETHER': [1],
              'HEAD': list(range(nid['PY_A'], len(order) + 1)), 'NPLUS': [nid['PY_A']], 'NITRILE': [nid['CN']],
              'VSITES': [s[0] for s in CHOL_VS]}
    for j, (r, bl) in enumerate(rings.items(), 1): groups[f'RING{j}'] = [nid[n] for n in bl]
    ring_names = {f'RING{j}': r for j, r in enumerate(rings, 1)}
    return dict(order=order, types=types, charges=charges, mass=mass, X=X, bonds=bonds, constr=constr,
                angles=angles, dihedrals=dihedrals, groups=groups, rings=rings, scale=scale, ring_names=ring_names,
                residue=M['residue'], nid=nid, beads=beads)

def write_itp(m, path):
    L = []
    L.append('; Martini 3 model of Chol4CSN, generated by cg_build.py')
    L.append('; Steroid: Martini 3 CHOL (Borges-Araujo et al. 2023), bead 1 (ROH) retyped as an ether bead.')
    L.append('; Linker/head: mapped from the all-atom mol2. Bead types and bonded values are STARTING values.')
    L.append('; Add  define = -DFLEXIBLE  in the minimisation mdp so constraints become stiff bonds.')
    L.append('')
    L.append('[ moleculetype ]\n; molname  nrexcl\n  %s     1\n' % m['residue'])
    L.append('[ atoms ]\n; id type resnr residue atom cgnr charge mass   ; AA atoms')
    mp = {b['name']: ' '.join(b['atoms']) for b in m['beads']}
    for i, (n, t, q, ms) in enumerate(zip(m['order'], m['types'], m['charges'], m['mass']), 1):
        mtxt = '' if ms is None else f'{ms:8.1f}'
        L.append(f'{i:3d} {t:<6s} 1 {m["residue"]} {n:<5s} {i:3d} {q:5.1f} {mtxt}   ; {mp.get(n, "CHOL template")}')
    L.append('\n[ bonds ]\n; i j funct length force')
    L.append('  %d %d 1 %.3f %d   ; CHOL C1-C2' % (CHOL_BOND[0], CHOL_BOND[1], CHOL_BOND[2], CHOL_BOND[3]))
    for b in m['bonds']:
        L.append(f'{b["i"]:3d} {b["j"]:3d} 1 {b["length"]:.3f} {b["k"]:g}   ; {m["order"][b["i"]-1]}-{m["order"][b["j"]-1]}')
    L.append('#ifdef FLEXIBLE')
    for (i, j, l) in CHOL_CONSTR: L.append(f'{i:3d} {j:3d} 1 {l:.5f} 100000')
    for c in m['constr']: L.append(f'{c["i"]:3d} {c["j"]:3d} 1 {c["length"]:.3f} 100000')
    L.append('#endif')
    L.append('\n[ constraints ]\n#ifndef FLEXIBLE')
    for (i, j, l) in CHOL_CONSTR: L.append(f'{i:3d} {j:3d} 1 {l:.5f}')
    for c in m['constr']: L.append(f'{c["i"]:3d} {c["j"]:3d} 1 {c["length"]:.3f}')
    L.append('#endif')
    L.append('\n[ angles ]\n; i j k funct theta force')
    L.append('  %d %d %d %d %.1f %d   ; CHOL R4-R3-C2' % CHOL_ANGLE)
    for a in m['angles']:
        L.append(f'{a["i"]:3d} {a["j"]:3d} {a["k"]:3d} 2 {a["theta"]:.1f} {a["k_force"]:g}   ; '
                 f'{m["order"][a["i"]-1]}-{m["order"][a["j"]-1]}-{m["order"][a["k"]-1]}')
    L.append('\n[ dihedrals ]\n; i j k l funct angle force')
    L.append('  %d %d %d %d %d %.1f %d   ; CHOL R6-R4-R3-C2' % CHOL_DIH)
    for d in m['dihedrals']:
        L.append(f'{d["i"]:3d} {d["j"]:3d} {d["k"]:3d} {d["l"]:3d} 2 {d["phi"]:.1f} {d.get("k_force", d.get("k")):g}'
                 + (f'   ; {d["note"]}' if d.get('note') else ''))
    L.append('\n[ virtual_sites3 ]')
    for s, a, b, c in CHOL_VS: L.append(f'{s:3d}   8  3  2  4   {a:.5f}   {b:.5f}   {c:.5f}')
    L.append('\n[ exclusions ]\n; steroid beads exclude each other (as in the CHOL model)')
    for i in range(1, 8): L.append(' '.join(str(x) for x in range(i, 9)))
    open(path, 'w').write('\n'.join(L) + '\n')

def write_gro(m, path):
    res = m['residue']; n = len(m['order'])
    with open(path, 'w') as f:
        f.write(f'{res} Martini 3\n{n}\n')
        for i, (nm, x) in enumerate(zip(m['order'], m['X']), 1):
            f.write(f'{1:5d}{res:<5s}{nm:>5s}{i:5d}{x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}\n')
        f.write('   10.00000   10.00000   10.00000\n')

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('mol2'); ap.add_argument('map'); ap.add_argument('--outdir', default='.')
    ap.add_argument('--overrides', default=None)
    a = ap.parse_args()
    M = json.load(open(a.map)); ov = json.load(open(a.overrides)) if a.overrides else None
    m = build(a.mol2, M, ov)
    os.makedirs(a.outdir, exist_ok=True)
    write_itp(m, os.path.join(a.outdir, 'CH4C.itp')); write_gro(m, os.path.join(a.outdir, 'CH4C_cg.gro'))
    json.dump(dict(natoms=len(m['order']), groups={k: dict(ids=v) for k, v in m['groups'].items()},
                   ring_names=m['ring_names']), open(os.path.join(a.outdir, 'cg_groups.json'), 'w'), indent=1)
    json.dump(dict(order=m['order'], bonds=m['bonds'], constr=m['constr'], angles=m['angles'], scale=m['scale'],
                   beads=m['beads'], ring_bond_nm=M['ring_bond_nm']), open(os.path.join(a.outdir, 'cg_index.json'), 'w'), indent=1)
    X = m['X']; n = len(m['order'])
    print(f'{n} beads, {len(m["bonds"])+1} bonds, {len(m["constr"])+3} constraints, {len(m["angles"])+1} angles')
    print('net charge %+.1f' % sum(m['charges']))
    ext = X.max(0) - X.min(0)
    print('bounding box (nm): %.2f %.2f %.2f' % tuple(ext))
    print('bonds (nm):', ' '.join(f'{m["order"][b["i"]-1]}-{m["order"][b["j"]-1]}={b["length"]:.3f}' for b in m['bonds']))
    ch = chol_template()
    D = np.linalg.norm(X[:, None] - X[None], axis=2); np.fill_diagonal(D, 9)
    bonded = {(b['i'], b['j']) for b in m['bonds']} | {(c['i'], c['j']) for c in m['constr']}
    mind = min(D[i, j] for i in range(n) for j in range(i + 1, n) if (i + 1, j + 1) not in bonded and i >= 8)
    print('closest non-bonded pair among head/linker beads and steroid: %.3f nm' % mind)
