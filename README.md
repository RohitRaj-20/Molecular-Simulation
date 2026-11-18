# Self-aggregation of Chol4CSN in water and DMSO (Martini 3)

**Rohit Raj** · Roll number 23ms176 · GitHub: [RohitRaj-20](https://github.com/RohitRaj-20)

## About

This project studies how a cholesterol-linked, positively charged α-cyanostilbene molecule (Chol4CSN, with bromide as the counter-ion) self-aggregates. Martini 3 coarse-grained simulations of the whole molecule are run in GROMACS, first in water and then in DMSO, with water/DMSO mixtures also included. Cluster sizes, ring stacking and bromide positions are analysed, and the potential of mean force (PMF) of two molecules approaching each other is calculated in each solvent to explain the differences in aggregation.

The molecule (130 atoms, charge +1) is mapped onto 24 Martini 3 beads. Nine come from the published Martini 3 cholesterol model and fifteen are mapped from the linker and head group. See `docs/runbook.md` for the full protocol and `docs/project_summary.md` for the short project description.

## Status

Set up and tested: model builder, topology checker, solvent-box, dimer and umbrella-window builders, aggregation analysis, PMF analysis. Not yet done: the GROMACS simulations. The bead types of the molecule are starting choices and still need testing (runbook Parts 5 and 9).

## Quick start

Install GROMACS 2021 or newer first (runbook Part 1.3), then:

```bash
git clone https://github.com/RohitRaj-20/chol4csn-martini.git ~/chol4csn_m3
cd ~/chol4csn_m3
python3 -m venv venv
source env.sh
pip install -r requirements.txt
bash ff/get_martini_ff.sh          # downloads the three Martini 3 files into ff/
bash tests/smoke_test.sh           # rebuilds the model and checks it
python3 scripts/check_setup.py     # checks GROMACS, Python modules and files
```

Then continue with `docs/runbook.md` from Part 3. Open a new terminal with `source ~/chol4csn_m3/env.sh`.

## Layout

| Folder | Contents |
|---|---|
| `docs/` | `runbook.md` (step-by-step protocol), `project_summary.md` |
| `params/` | cleaned molecule file, bead mapping (`chol4csn_cg_map.json`), generated `CH4C.itp` and `CH4C_cg.gro` |
| `mdp/` | Martini 3 settings: `em`, `eq1`, `eq2`, `md_prod`, `pull_steered` |
| `scripts/` | model builder, topology checker, solvent/dimer/index/window builders, analysis scripts |
| `ff/` | Martini 3 download script (the force-field files themselves are not stored) |
| `tests/` | `smoke_test.sh` |
| `results/` | small final results (tables, figures) |

Simulation folders (`solvents/`, `dimer/`, `agg/`, `swap/`, `test/`) are created while running and are ignored by git, together with trajectories and run files.

## Workflow

1. Build the Martini model (`scripts/cg_build.py`), check it (`scripts/check_cg_itp.py`).
2. Build and equilibrate water, DMSO and three water/DMSO mixtures (`scripts/make_solvent_cg.py`).
3. Validate one molecule alone in water and in DMSO.
4. Dimer PMF per solvent: steered pull, umbrella windows (`scripts/setup_windows_cg.py`), `gmx wham`, `scripts/pmf_analysis.py`.
5. Aggregation: 24 molecules and 24 Br⁻ in a 10 nm box, three replicates per solvent, then solvent-swap tests.
6. Analysis: `scripts/analyze_agg_cg.py`, `gmx rdf`, summaries and plots.

## Limits

Martini 3 uses one dielectric constant (ε_r = 15) in every solvent, so solvent differences come only from the bead types, and bromide behaviour in DMSO is less reliable. The PMF is a coarse-grained free energy and is best used to compare solvents. See the end of the runbook for the full list.

## Citing

Please cite the Martini 3 force field (Souza et al., Nature Methods, 2021), the Martini 3 cholesterol model (Borges-Araújo et al., 2023) and GROMACS. See `NOTICE.md`.

## Licence

MIT for the code and documents in this repository (see `LICENSE`). Third-party material is listed in `NOTICE.md`.
