#!/usr/bin/env python3
"""Check that everything the runbook needs is in place. Run it after Part 2 and again before every new stage.
Usage:  python3 check_setup.py            (uses $ROOT, or ~/chol4csn_m3)
Prints one line per check: OK / MISSING / WARN, and what to do about it."""
import os, shutil, subprocess, importlib, sys

ROOT = os.environ.get('ROOT', os.path.expanduser('~/chol4csn_m3'))
bad = 0
def report(kind, msg, fix=''):
    global bad
    if kind == 'MISSING': bad += 1
    print(f'{kind:8s}{msg}' + (f'\n         -> {fix}' if fix and kind != 'OK' else ''))

gmx = shutil.which('gmx') or shutil.which('gmx_mpi')
if gmx:
    out = subprocess.run([gmx, '--version'], capture_output=True, text=True).stdout
    ver = next((l.split(':')[1].strip() for l in out.splitlines() if l.startswith('GROMACS version')), '?')
    gpu = next((l.split(':')[1].strip() for l in out.splitlines() if l.startswith('GPU support')), '?')
    report('OK', f'GROMACS {ver} at {gmx}; GPU support: {gpu}')
    try:
        major = int(ver.split('.')[0].replace('VERSION', '').strip()[:4])
        if major < 2021: report('WARN', 'GROMACS older than 2021', 'install a newer release (runbook Part 1.3)')
    except ValueError: pass
    if gpu.lower() in ('disabled', 'none', '?'):
        report('WARN', 'this GROMACS has no GPU support', 'fine for the tests; remove "-nb gpu" from MDRUN_OPTS in env.sh, or build a GPU version for production')
else:
    report('MISSING', 'gmx not found on PATH', 'install GROMACS (runbook Part 1.3) and run: source ~/chol4csn_m3/env.sh')

for m in ('numpy', 'scipy', 'matplotlib'):
    try: importlib.import_module(m); report('OK', f'python module {m}')
    except ImportError: report('MISSING', f'python module {m}', f'pip install {m}   (inside the venv)')
try: importlib.import_module('MDAnalysis'); report('OK', 'python module MDAnalysis (optional, only for .xtc analysis)')
except ImportError: report('WARN', 'MDAnalysis not installed (optional)', 'pip install MDAnalysis, or analyse .gro trajectories instead')

need = {
 'ff/martini_v3.0.0.itp': 'Part 2.2', 'ff/martini_v3.0.0_ions_v1.itp': 'Part 2.2', 'ff/martini_v3.0.0_solvents_v1.itp': 'Part 2.2',
 'params/Chol4CSN_clean.mol2': 'Part 2.1', 'params/chol4csn_cg_map.json': 'Part 2.1',
 'scripts/cg_build.py': 'Part 2.1', 'scripts/make_solvent_cg.py': 'Part 2.1', 'scripts/setup_windows_cg.py': 'Part 2.1',
 'mdp/em.mdp': 'Part 2.1', 'mdp/eq1.mdp': 'Part 2.1', 'mdp/eq2.mdp': 'Part 2.1', 'mdp/md_prod.mdp': 'Part 2.1', 'mdp/pull_steered.mdp': 'Part 2.1',
}
for f, part in need.items():
    p = os.path.join(ROOT, f)
    report('OK' if os.path.exists(p) else 'MISSING', f, f'see runbook {part}; expected at {p}')
print('\nAll required items present.' if not bad else f'\n{bad} required item(s) missing: fix them, then run this again.')
sys.exit(1 if bad else 0)
