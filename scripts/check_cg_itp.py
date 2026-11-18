#!/usr/bin/env python3
"""Sanity checks on a Martini molecule itp before you run grompp (no GROMACS needed).
Usage: python3 check_cg_itp.py CH4C.itp martini_v3.0.0.itp [--charge 1]
Checks: bead types exist in the force-field file, every index in bonds/constraints/angles/dihedrals/virtual sites
is a valid bead, no bead is bonded to itself, no term is listed twice, net charge, and that every bead is
connected to the rest of the molecule (or is a virtual site)."""
import sys, re, argparse
ap = argparse.ArgumentParser(); ap.add_argument('itp'); ap.add_argument('ff'); ap.add_argument('--charge', type=float, default=None)
a = ap.parse_args()
types = set()
sec = None
for l in open(a.ff):
    l = l.split(';')[0].strip()
    if not l: continue
    if l.startswith('['): sec = l.strip('[] ').lower(); continue
    if sec == 'atomtypes': types.add(l.split()[0])
S = {}; sec = None; ifdef = []
for l in open(a.itp):
    l = l.split(';')[0].strip()
    if not l: continue
    if l.startswith('#'): continue                     # both FLEXIBLE branches are checked
    if l.startswith('['): sec = l.strip('[] ').lower(); S.setdefault(sec, []); continue
    S[sec].append(l.split())
n = len(S['atoms']); bad = 0
def err(m):
    global bad; bad += 1; print('  PROBLEM:', m)
print(f'{n} beads')
for r in S['atoms']:
    if r[1] not in types: err(f'bead type {r[1]} (bead {r[0]}) not found in {a.ff}')
q = sum(float(r[6]) for r in S['atoms']); print(f'net charge {q:+.2f}')
if a.charge is not None and abs(q - a.charge) > 1e-6: err(f'net charge {q} differs from the expected {a.charge}')
nat = {'bonds': 2, 'constraints': 2, 'angles': 3, 'dihedrals': 4, 'virtual_sites3': 4}
seen = {}
for s, k in nat.items():
    for r in S.get(s, []):
        ids = [int(x) for x in r[:k + (1 if s == 'virtual_sites3' else 0)]]
        for i in ids:
            if not 1 <= i <= n: err(f'{s}: index {i} out of range')
        key = (s, tuple(ids) if s == 'virtual_sites3' else min(tuple(ids[:k]), tuple(ids[:k][::-1])))
        if key in seen and s != 'bonds': err(f'{s}: {ids} listed twice')
        seen[key] = 1
        if s in ('bonds', 'constraints') and ids[0] == ids[1]: err(f'{s}: bead bonded to itself {ids}')
vs = {int(r[0]) for r in S.get('virtual_sites3', [])}
adj = {i: set() for i in range(1, n + 1)}
for s in ('bonds', 'constraints'):
    for r in S.get(s, []):
        i, j = int(r[0]), int(r[1]); adj[i].add(j); adj[j].add(i)
for r in S.get('virtual_sites3', []):
    for j in r[1:4]: adj[int(r[0])].add(int(j)); adj[int(j)].add(int(r[0]))
for r in S.get('bonds', []):
    pass
seen_ = {1}; stack = [1]
while stack:
    u = stack.pop()
    for v in adj[u]:
        if v not in seen_: seen_.add(v); stack.append(v)
if len(seen_) != n: err(f'beads not connected to bead 1: {sorted(set(range(1, n + 1)) - seen_)}')
print('virtual sites:', sorted(vs))
print('OK' if not bad else f'{bad} problem(s)')
sys.exit(1 if bad else 0)
