#!/usr/bin/env python3
"""Martini 3 version of setup_windows.py: dt = 20 fs, reaction-field electrostatics, Parrinello-Rahman tau_p = 12.
Pick start frames from the steered pull and write one umbrella window per target distance.
Usage:
  python3 setup_windows.py --pullx pull_pullx.xvg --xtc pull.xtc --tpr pull.tpr \
      --ndx index.ndx --top topol.top --gro-ref dimer.gro --outdir windows [--dry-run]
Creates windows/win_000/ ... each with conf.gro, umbrella.mdp, and (unless --dry-run) umb.tpr.
Also writes windows/windows.csv, tpr-files.dat and pullx-files.dat for gmx wham.
Segments are start:stop:step:k  (nm, nm, nm, kJ/mol/nm^2)."""
import argparse, os, subprocess, csv
import numpy as np

SEGMENTS = "0.5:1.5:0.05:1000,1.6:3.0:0.1:500,3.2:4.4:0.2:250"

MDP = """; Martini 3 umbrella window  r0 = {r0:.3f} nm   k = {k:g}
integrator               = md
dt                       = 0.02
nsteps                   = {nsteps}          ; {ns:g} ns
nstlog                   = 5000
nstenergy                = 5000
nstxout-compressed       = 1000              ; every 20 ps
compressed-x-grps        = SOLUTE_BR         ; solutes and ions only, keeps files small
gen-vel                  = yes
gen-temp                 = 300
gen-seed                 = -1
continuation             = no
constraints              = none
constraint-algorithm     = LINCS
lincs-order              = 12
lincs-iter               = 2
cutoff-scheme            = Verlet
verlet-buffer-tolerance  = 0.005
coulombtype              = reaction-field
rcoulomb                 = 1.1
epsilon-r                = 15
epsilon-rf               = 0
vdwtype                  = cutoff
vdw-modifier             = potential-shift-verlet
rvdw                     = 1.1
tcoupl                   = v-rescale
tc-grps                  = System
tau-t                    = 1.0
ref-t                    = 300
pcoupl                   = parrinello-rahman
pcoupltype               = isotropic
tau-p                    = 12.0
ref-p                    = 1.0
compressibility          = 3e-4
pbc                      = xyz
pull                     = yes
pull-ngroups             = 2
pull-ncoords             = 1
pull-group1-name         = MOL1
pull-group2-name         = MOL2
pull-group1-pbcatom      = {pbc1}
pull-group2-pbcatom      = {pbc2}
pull-coord1-type         = umbrella
pull-coord1-geometry     = distance
pull-coord1-groups       = 1 2
pull-coord1-dim          = Y Y Y
pull-coord1-start        = no
pull-coord1-init         = {r0:.3f}
pull-coord1-rate         = 0.0
pull-coord1-k            = {k:g}
pull-nstxout             = 50                ; every 1 ps
pull-nstfout             = 50
"""

def parse_segments(spec):
    out = []
    for seg in spec.split(','):
        a, b, s, k = (float(x) for x in seg.split(':'))
        n = int(round((b - a) / s)) + 1
        out += [(round(a + i * s, 4), k) for i in range(n)]
    return out

def read_xvg(path):
    rows = [l.split() for l in open(path) if l.strip() and l[0] not in '#@']
    arr = np.array([[float(x) for x in r] for r in rows])
    return arr[:, 0], arr[:, 1]          # time (ps), pull coordinate (nm)

def pick_frames(t, d, targets):
    """For every target distance return (time, actual distance) of the nearest frame."""
    return [(float(t[np.argmin(np.abs(d - r))]), float(d[np.argmin(np.abs(d - r))])) for r in targets]

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--pullx', required=True); ap.add_argument('--xtc'); ap.add_argument('--tpr')
    ap.add_argument('--ndx'); ap.add_argument('--top'); ap.add_argument('--outdir', default='windows')
    ap.add_argument('--segments', default=SEGMENTS)
    ap.add_argument('--ns', type=float, default=500.0, help='length of each window, ns (Martini: 0.5 us is a good first pass)')
    ap.add_argument('--pbc1', type=int, required=True); ap.add_argument('--pbc2', type=int, required=True)
    ap.add_argument('--dry-run', action='store_true')
    a = ap.parse_args()

    t, d = read_xvg(a.pullx)
    wins = parse_segments(a.segments)
    frames = pick_frames(t, d, [w[0] for w in wins])
    os.makedirs(a.outdir, exist_ok=True)
    rows, tprs, pxs = [], [], []
    for i, ((r0, k), (tt, dd)) in enumerate(zip(wins, frames)):
        wd = os.path.join(a.outdir, f'win_{i:03d}'); os.makedirs(wd, exist_ok=True)
        nsteps = int(a.ns * 1000 / 0.02)
        open(os.path.join(wd, 'umbrella.mdp'), 'w').write(
            MDP.format(r0=r0, k=k, nsteps=nsteps, ns=a.ns, pbc1=a.pbc1, pbc2=a.pbc2))
        flag = ''
        if abs(dd - r0) > 0.03: flag = 'FRAME FAR FROM TARGET'
        rows.append([i, r0, k, f'{tt:.1f}', f'{dd:.3f}', flag])
        tprs.append(f'win_{i:03d}/umb.tpr'); pxs.append(f'win_{i:03d}/umb_pullx.xvg')
        if not a.dry_run:
            subprocess.run(['gmx', 'trjconv', '-s', a.tpr, '-f', a.xtc, '-dump', f'{tt}',
                            '-o', os.path.join(wd, 'conf.gro')], input='0\n', text=True, check=True)
            subprocess.run(['gmx', 'grompp', '-f', os.path.join(wd, 'umbrella.mdp'),
                            '-c', os.path.join(wd, 'conf.gro'), '-p', a.top, '-n', a.ndx,
                            '-o', os.path.join(wd, 'umb.tpr'), '-maxwarn', '0'], check=True)
    with open(os.path.join(a.outdir, 'windows.csv'), 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['window', 'r0_nm', 'k', 'time_ps', 'frame_dist_nm', 'note']); w.writerows(rows)
    open(os.path.join(a.outdir, 'tpr-files.dat'), 'w').write('\n'.join(tprs) + '\n')
    open(os.path.join(a.outdir, 'pullx-files.dat'), 'w').write('\n'.join(pxs) + '\n')
    print(f'{len(wins)} windows written to {a.outdir}')
    bad = [r for r in rows if r[5]]
    if bad: print('check these windows (no steered frame near the target):', [r[0] for r in bad])
