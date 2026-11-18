#!/usr/bin/env python3
"""Mean +/- SEM across replicate CSVs from analyze_agg.py.
Usage: python3 summarize_reps.py agg_rep1.csv agg_rep2.csv agg_rep3.csv [--last-ns 100] [--nmol 24]"""
import sys, argparse, csv
import numpy as np
ap = argparse.ArgumentParser(); ap.add_argument('csvs', nargs='+')
ap.add_argument('--last-ns', type=float, default=100.0); ap.add_argument('--nmol', type=int, default=24)
a = ap.parse_args()
cols = ['n_clusters', 'n_monomers', 'largest', 'mean_size_w', 'rg_largest_nm', 'kappa2_largest',
        'stacked_pairs', 'mol_in_stacks', 'br_per_cation']
res = {c: [] for c in cols}; t_dimer = []; t_half = []
for p in a.csvs:
    rows = list(csv.DictReader(open(p)))
    t = np.array([float(r['time_ps']) for r in rows]) / 1000.0         # ns
    def col(c): return np.array([float(r[c]) if r[c] != 'nan' else np.nan for r in rows])
    m = t >= t.max() - a.last_ns
    for c in cols: res[c].append(np.nanmean(col(c)[m]))
    big = col('largest')
    t_dimer.append(t[big >= 2][0] if (big >= 2).any() else np.nan)
    t_half.append(t[big >= a.nmol / 2][0] if (big >= a.nmol / 2).any() else np.nan)
def ms(v):
    v = np.array(v, float); n = np.sum(~np.isnan(v))
    return np.nanmean(v), (np.nanstd(v, ddof=1) / np.sqrt(n) if n > 1 else np.nan), int(n)
print(f'{len(a.csvs)} replicates, averages over the last {a.last_ns:g} ns of each')
for c in cols: m, s, n = ms(res[c]); print(f'{c:16s} {m:8.3f} +/- {s:6.3f}')
for name, v in (('t_first_dimer_ns', t_dimer), ('t_half_in_one_cluster_ns', t_half)):
    m, s, n = ms(v); print(f'{name:24s} {m:8.1f} +/- {s:6.1f}   ({n} of {len(v)} replicates reached it)')
