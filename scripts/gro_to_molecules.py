#!/usr/bin/env python3
"""Print the [ molecules ] block for a .gro file, in file order.
Usage: python3 gro_to_molecules.py system.gro [RES=MOLNAME ...]
Molecules are counted by changes of residue number/name, so every molecule needs its own residue number.
Martini example (residue names in the gro differ from the molecule names in the topology):
  python3 gro_to_molecules.py system.gro ION=BR DMS=DMSO
Multi-bead solvents (DMSO) are counted once per residue, as they should be."""
import sys
args = [x for x in sys.argv[1:] if '=' not in x]
rename = dict(x.split('=') for x in sys.argv[1:] if '=' in x)
lines = open(args[0]).read().splitlines()
n = int(lines[1]); prev = None; blocks = []
for l in lines[2:2 + n]:
    key = (l[0:5], l[5:10].strip())
    if key != prev:
        nm = rename.get(key[1], key[1])
        if blocks and blocks[-1][0] == nm: blocks[-1][1] += 1
        else: blocks.append([nm, 1])
        prev = key
print('[ molecules ]')
for name, c in blocks: print(f'{name:<10s} {c}')
