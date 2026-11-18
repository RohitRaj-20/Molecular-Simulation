#!/usr/bin/env python3
"""Quick plots.
  python3 plot_results.py pmf  --histo windows/histo.xvg --pmf pmf_physical.csv [--pmf2 other.csv ...] -o pmf.png
  python3 plot_results.py agg  agg_rep1.csv agg_rep2.csv agg_rep3.csv -o agg.png
'pmf' shows the umbrella histograms (overlap check) above the physical PMF(s).
'agg'  shows largest cluster, monomers and ring-stacked pairs against time, one line per replicate."""
import sys, argparse, csv
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

def read_xvg(p):
    rows = [l.split() for l in open(p) if l.strip() and l[0] not in '#@']
    return np.array([[float(x) for x in r] for r in rows])

ap = argparse.ArgumentParser(); ap.add_argument('mode', choices=['pmf', 'agg'])
ap.add_argument('csvs', nargs='*'); ap.add_argument('--histo'); ap.add_argument('--pmf'); ap.add_argument('--pmf2', nargs='*', default=[])
ap.add_argument('-o', default='plot.png'); a = ap.parse_args()
if a.mode == 'pmf':
    fig, ax = plt.subplots(2, 1, figsize=(7, 7), sharex=True)
    h = read_xvg(a.histo)
    for j in range(1, h.shape[1]): ax[0].plot(h[:, 0], h[:, j], lw=0.8)
    ax[0].set_ylabel('umbrella histograms (counts)')
    for p in [a.pmf] + a.pmf2:
        d = np.loadtxt(p, delimiter=',', skiprows=1); ax[1].plot(d[:, 0], d[:, 1], label=p)
    ax[1].set_xlabel('COM distance (nm)'); ax[1].set_ylabel('W(r) (kJ/mol)'); ax[1].axhline(0, c='k', lw=0.5); ax[1].legend()
else:
    fig, ax = plt.subplots(3, 1, figsize=(7, 8), sharex=True)
    for p in a.csvs:
        rows = list(csv.DictReader(open(p))); t = np.array([float(r['time_ps']) for r in rows]) / 1000
        for k, c in enumerate(['largest', 'n_monomers', 'stacked_pairs']):
            ax[k].plot(t, [float(r[c]) for r in rows], lw=0.9, label=p)
    for k, c in enumerate(['largest cluster (molecules)', 'monomers', 'stacked ring pairs']): ax[k].set_ylabel(c)
    ax[2].set_xlabel('time (ns)'); ax[0].legend(fontsize=7)
fig.tight_layout(); fig.savefig(a.o, dpi=150); print('wrote', a.o)
