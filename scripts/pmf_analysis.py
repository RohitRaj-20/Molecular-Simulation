#!/usr/bin/env python3
"""Turn gmx wham output into a physical PMF and a standard-state binding free energy.
Usage:
  python3 pmf_analysis.py profile.xvg [--bsprof bsProfs.xvg] [--plateau 3.8 4.4] [--rb 1.0] [--T 300]
Steps
  1. W(r) = W_wham(r) + 2 kT ln r          (removes the 4 pi r^2 volume term)
  2. shift so the mean of W over the plateau window is zero
  3. K = 4 pi * integral_0^rb r^2 exp(-W/kT) dr   (nm^3), dG0 = -kT ln(K * C0), C0 = 1 M = 0.6022 nm^-3
  4. if bootstrap profiles are given, repeat 1-3 on each and report the spread
rb is the outer edge of the bound state; try several (e.g. the barrier position)."""
import sys, argparse
import numpy as np

KB = 0.0083144626   # kJ/mol/K
C0 = 0.602214       # molecules per nm^3 for 1 mol/L

def read_xvg(path):
    rows = [l.split() for l in open(path) if l.strip() and l[0] not in '#@']
    return np.array([[float(x) for x in r] for r in rows])

def physical_pmf(r, w, T, plateau):
    kT = KB * T
    W = w + 2 * kT * np.log(r)
    m = (r >= plateau[0]) & (r <= plateau[1])
    if not m.any(): raise SystemExit('no points in the plateau window; widen --plateau')
    return W - W[m].mean()

def dG0(r, W, T, rb, r_lo=0.0):
    kT = KB * T
    m = (r <= rb)
    rr, ww = r[m], W[m]
    # extend the first bin down to r_lo is not needed: the integrand is ~0 below the steric wall
    _trap = getattr(np, 'trapezoid', None) or np.trapz
    K = 4 * np.pi * _trap(rr ** 2 * np.exp(-ww / kT), rr)
    return -kT * np.log(K * C0)

def summarize(r, W, T, rb):
    i = int(np.argmin(W[r <= rb])) if (r <= rb).any() else int(np.argmin(W))
    return dict(Wmin=W[i], rmin=r[i], dG0=dG0(r, W, T, rb))

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('profile'); ap.add_argument('--bsprof')
    ap.add_argument('--plateau', nargs=2, type=float, default=[3.8, 4.4])
    ap.add_argument('--rb', type=float, nargs='+', default=[1.0, 1.5, 2.0])
    ap.add_argument('--T', type=float, default=300.0)
    ap.add_argument('--out', default='pmf_physical.csv')
    a = ap.parse_args()
    d = read_xvg(a.profile); r, w = d[:, 0], d[:, 1]
    if r.max() > 20: sys.exit('x values look like Angstrom or ps; expected nm')
    W = physical_pmf(r, w, a.T, a.plateau)
    np.savetxt(a.out, np.c_[r, W], delimiter=',', header='r_nm,W_kJ_per_mol', comments='')
    bs = None
    if a.bsprof:
        b = read_xvg(a.bsprof); assert np.allclose(b[:, 0], r), 'bootstrap x grid differs from profile'
        bs = [physical_pmf(r, b[:, j], a.T, a.plateau) for j in range(1, b.shape[1])]
    print(f'{"rb/nm":>6} {"Wmin":>9} {"rmin/nm":>8} {"dG0":>9}' + ('   +/- (bootstrap)' if bs else ''))
    for rb in a.rb:
        s = summarize(r, W, a.T, rb); line = f'{rb:6.2f} {s["Wmin"]:9.2f} {s["rmin"]:8.3f} {s["dG0"]:9.2f}'
        if bs:
            ss = [summarize(r, x, a.T, rb) for x in bs]
            line += f'   {np.std([q["dG0"] for q in ss]):.2f}  (Wmin +/- {np.std([q["Wmin"] for q in ss]):.2f})'
        print(line)
    print('units kJ/mol; wrote', a.out)
