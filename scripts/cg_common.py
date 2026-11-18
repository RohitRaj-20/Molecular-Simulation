"""Small helpers shared by the Martini scripts (numpy only)."""
import numpy as np

DEFAULT_MASS = {'T': 36.0, 'S': 54.0}      # Martini 3 default masses by bead size (regular = 72)

def parse_itp_atoms(path):
    """Return names, types, charges, masses of the first [ atoms ] block of a Martini itp."""
    names, types, charges, masses = [], [], [], []
    inatoms = False
    for line in open(path):
        l = line.split(';')[0].strip()
        if not l: continue
        if l.startswith('['):
            inatoms = l.replace(' ', '').lower().startswith('[atoms]')
            if names and not inatoms: break
            continue
        if inatoms:
            f = l.split()
            names.append(f[4]); types.append(f[1]); charges.append(float(f[6]))
            masses.append(float(f[7]) if len(f) > 7 else DEFAULT_MASS.get(f[1][0], 72.0))
    return names, types, np.array(charges), np.array(masses)

def read_gro(path):
    """Single-frame gro: names, residue names, positions (nm), box (3,)."""
    L = open(path).read().splitlines(); n = int(L[1])
    resn = [l[5:10].strip() for l in L[2:2 + n]]; names = [l[10:15].strip() for l in L[2:2 + n]]
    xyz = np.array([[float(l[20:28]), float(l[28:36]), float(l[36:44])] for l in L[2:2 + n]])
    box = np.array([float(v) for v in L[2 + n].split()[:3]])
    return names, resn, xyz, box

def gro_frames(path):
    """Generator over frames of a multi-frame gro: (time_ps or frame index, xyz (nm), box (3,))."""
    with open(path) as f:
        k = 0
        while True:
            title = f.readline()
            if not title: return
            n = int(f.readline())
            rows = [f.readline() for _ in range(n)]
            box = np.array([float(v) for v in f.readline().split()[:3]])
            xyz = np.array([[float(l[20:28]), float(l[28:36]), float(l[36:44])] for l in rows])
            t = float(title.split('t=')[1].split()[0]) if 't=' in title else float(k)
            yield t, xyz, box
            k += 1
