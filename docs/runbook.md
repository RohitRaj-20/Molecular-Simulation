> **Repository note.** This repository is the project folder the runbook calls `~/chol4csn_m3` (`$ROOT`). Clone it into `~/chol4csn_m3` as described in the repository README, and then: skip Part 1.2 (`env.sh` is already included), skip Part 2.1 (the files are already here) and, in Part 2.2, you can run `bash ff/get_martini_ff.sh` instead of the manual download. Everything else in the runbook applies unchanged.

# Chol4CSN in Martini 3: complete start-to-finish runbook

This version assumes you have done nothing yet. It starts from an empty computer and ends with the dimer PMFs and the aggregation simulations in water and DMSO, all with the Martini 3 coarse-grained force field. Martini is a force field, not a program, so it runs inside a simulation engine; this runbook uses GROMACS, the standard engine for Martini 3.

Your molecule: a cholesteryl ether with a four-carbon linker to a cationic head group (pyridinium–phenyl–α-cyanostilbene–phenyl–pyridine), 130 atoms, charge +1, with Br⁻ as the counter-ion. The whole molecule is simulated. It becomes 24 Martini beads.

---

## Part 0. Read this first

### 0.1 What you will do, in order

| Stage | What | Result | Gate (do not go on until true) |
|---|---|---|---|
| 1 | Prepare the computer: Python, GROMACS | `gmx` and `python3` work | `check_setup.py` shows no MISSING |
| 2 | Get the project files and the Martini 3 force field | folder `~/chol4csn_m3` filled | `check_setup.py` shows no MISSING |
| 3 | Build the Martini model of your molecule | `CH4C.itp`, `CH4C_cg.gro` | `check_cg_itp.py` prints OK |
| 4 | Build and equilibrate the solvent boxes | water, DMSO, 3 mixtures | densities as expected |
| 5 | Run one molecule alone in water and DMSO | molecule is stable | passes the checks in 5.5 |
| 6 | Dimer PMF in each solvent | free-energy curves | window histograms overlap; halves agree |
| 7 | Aggregation simulations (24 molecules) | trajectories | runs finish |
| 8 | Analyse | cluster sizes, RDFs, plots | replicates agree |
| 9 | Controls and later improvements | confidence in the answer | |

Start with water only. Do stages 3–8 for water completely, including a short first-pass PMF, before touching DMSO.

### 0.2 What you need

- A Linux computer (Ubuntu 22.04 or newer is the easiest), or a Mac, or Windows with WSL2. A university cluster also works.
- A recent NVIDIA GPU is strongly recommended for stages 6 and 7. Stages 1–5 run on a CPU.
- About 100 GB of free disk for all solvents (more if you keep every trajectory).
- The two files from me: `chol4csn_martini_scripts.zip` and this document.

### 0.3 How to read the commands

- Text in grey boxes is typed into a terminal, line by line, in order. Lines starting with `#` are comments; you do not type them.
- `$ROOT`, `$S`, `$MDP` and `$MDRUN_OPTS` are shortcuts defined once in Part 1.2. Every time you open a new terminal, run `source ~/chol4csn_m3/env.sh` first.
- After each command block, "Expect" says what correct output looks like. If yours differs, stop and read Troubleshooting at the end.
- Never add `-maxwarn` to `gmx grompp`. A warning means something is wrong.
- Long runs: start them inside `tmux` (`tmux new -s md`; detach with Ctrl-b then d; return with `tmux attach -t md`) so they survive closing the terminal.

### 0.4 Honest status of this material

I could not run GROMACS here, so every `gmx` command is untested. The Python scripts were tested on your mol2 and on synthetic data with known answers, and the generated topology passed `check_cg_itp.py` against the real `martini_v3.0.0.itp`. The molecule's bead types are my starting choices, not validated parameters; Parts 5 and 9 tell you how to test them.

### 0.5 Words you will meet

| Word | Meaning |
|---|---|
| Bead | one coarse-grained particle; a regular Martini bead stands for about 4 heavy atoms, a tiny bead for about 2 |
| W | the Martini water bead; one W stands for 4 water molecules |
| `.gro` | coordinate file (positions in nm) |
| `.itp` / `.top` | topology: which beads exist and how they are connected |
| `.mdp` | the settings file for one simulation step |
| `grompp` | builds a run input file `.tpr` from `.mdp` + coordinates + topology |
| `mdrun` | runs the simulation |
| PMF | potential of mean force: free energy as a function of the distance between the two molecules |
| Umbrella sampling | many simulations ("windows"), each holding the two molecules at a different distance with a spring |
| WHAM | the method that joins all windows into one PMF (`gmx wham`) |
| Replicate | an independent repeat with a different random seed |

---

## Part 1. Prepare the computer

### 1.1 Folders and Python

```bash
mkdir -p ~/chol4csn_m3 && cd ~/chol4csn_m3
sudo apt update && sudo apt install -y python3 python3-venv python3-pip git unzip tmux      # Ubuntu / WSL2
# Mac: install Homebrew (brew.sh), then:  brew install python git
python3 -m venv ~/chol4csn_m3/venv
source ~/chol4csn_m3/venv/bin/activate
pip install numpy scipy matplotlib
python3 -c "import numpy, scipy, matplotlib; print('python ok')"
```
Expect: `python ok`. On a cluster you may not have sudo; use the cluster's Python module instead (`module avail python`) and skip the `apt` line.

### 1.2 The environment file

```bash
cat > ~/chol4csn_m3/env.sh <<'EOT'
export ROOT=$HOME/chol4csn_m3
export S=$ROOT/scripts
export MDP=$ROOT/mdp
source $ROOT/venv/bin/activate
[ -f $HOME/gromacs/bin/GMXRC ] && source $HOME/gromacs/bin/GMXRC
# mdrun options: edit -ntomp to your CPU cores (type: nproc). Remove "-nb gpu" if you have no GPU.
export MDRUN_OPTS="-ntmpi 1 -ntomp 8 -nb gpu"
EOT
source ~/chol4csn_m3/env.sh
echo $ROOT
```
Expect: `/home/<you>/chol4csn_m3`. Use `source ~/chol4csn_m3/env.sh` at the start of every terminal session from now on.

### 1.3 Install GROMACS (version 2021 or newer; newest is best)

Pick one route.

**A. University cluster.** GROMACS is usually installed already:
```bash
module avail gromacs
module load gromacs/<version>        # put this line in env.sh as well
gmx --version | head -n 12
```

**B. Ubuntu or WSL2 workstation with an NVIDIA GPU (recommended for production).**
Check the GPU and CUDA first:
```bash
nvidia-smi                           # must list your GPU
nvcc --version                       # CUDA toolkit; if missing: sudo apt install nvidia-cuda-toolkit
sudo apt install -y build-essential cmake
```
Download the newest release from gromacs.org (Downloads page) and replace `VERSION` below by its number:
```bash
cd ~ && wget https://ftp.gromacs.org/gromacs/gromacs-VERSION.tar.gz
tar xfz gromacs-VERSION.tar.gz && cd gromacs-VERSION && mkdir build && cd build
cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=OFF -DGMX_GPU=CUDA -DCMAKE_INSTALL_PREFIX=$HOME/gromacs
make -j 8
make install
source $HOME/gromacs/bin/GMXRC
gmx --version | head -n 12
```
Expect in the output: `GROMACS version: ...` and `GPU support: CUDA`. The build takes 10–30 minutes. `env.sh` already sources `~/gromacs/bin/GMXRC`.

**C. No GPU, or a quick start (CPU only).**
```bash
# Ubuntu:  sudo apt install -y gromacs
# Mac:     brew install gromacs
# any OS with conda:  conda install -c conda-forge gromacs
gmx --version | head -n 6
```
A CPU build is enough for Parts 3–5 and for testing everything. For production, move to route A or B. In `env.sh` change `MDRUN_OPTS` to `"-ntmpi 1 -ntomp 8"` (no `-nb gpu`).

If you only have a package-manager version, check its number (`gmx --version`); older than 2021 is too old.

### 1.4 A viewer

Install VMD (ks.uiuc.edu/Research/vmd, free) or PyMOL, to look at `.gro` files. You will use it in 3.4.

---

## Part 2. Get the project files

### 2.1 Unpack my zip

Move `chol4csn_martini_scripts.zip` to `~` (your home folder), then:
```bash
cd ~
unzip -q chol4csn_martini_scripts.zip
cp -r ~/chol4csn_martini/* ~/chol4csn_m3/ && rm -r ~/chol4csn_martini
mkdir -p ~/chol4csn_m3/{ff,aa_ref,solvents,test,dimer,agg,swap}
source ~/chol4csn_m3/env.sh
ls $ROOT
```
Expect: `agg  aa_ref  dimer  env.sh  ff  mdp  params  scripts  solvents  swap  test  venv` (plus `README.txt`).

The zip contains the scripts, the five mdp files, your molecule file (cleaned copy `Chol4CSN_clean.mol2`: residue name CH4C, atom names C1–C56, H1–H70, N1–N3, O1, charge +1) and the bead mapping (`chol4csn_cg_map.json`). Always use the cleaned mol2, not the original, because the mapping refers to its atom names.

### 2.2 Download the Martini 3 force field

```bash
cd $ROOT/ff
git clone --depth 1 https://github.com/Martini-Force-Field-Initiative/M3-Lipid-Parameters.git
cp M3-Lipid-Parameters/ITPs/martini_v3.0.0.itp \
   M3-Lipid-Parameters/ITPs/martini_v3.0.0_ions_v1.itp \
   M3-Lipid-Parameters/ITPs/martini_v3.0.0_solvents_v1.itp .
ls *.itp
grep -n -A4 "BR" martini_v3.0.0_ions_v1.itp | head -n 8
grep -n -A6 "DMSO" martini_v3.0.0_solvents_v1.itp | head -n 10
```
Expect: three `.itp` files; `BR` with bead type `SQ4` and charge −1; `DMSO` with two beads `S1` (SC6) and `O2` (TP6). If `git clone` fails, download the same three files from the Martini 3 downloads page at cgmartini.nl.

Cite the Martini 3 force-field paper (Souza et al., Nat. Methods 2021) and, because your cholesterol part reuses it, the Martini 3 cholesterol paper (Borges-Araújo et al., 2023).

### 2.3 Check everything is in place

```bash
python3 $S/check_setup.py
```
Expect: every line `OK` and the last line `All required items present.` (`MDAnalysis` may show `WARN`; that is optional.) Fix anything marked `MISSING` before going on.

---

## Part 3. Build the Martini model of your molecule

### 3.1 Generate the model

```bash
cd $ROOT/params
python3 $S/cg_build.py Chol4CSN_clean.mol2 chol4csn_cg_map.json --outdir .
```
Expect: `24 beads, 8 bonds, 15 constraints, 19 angles`, `net charge +1.0`, and a bounding box with a long side of about 3 nm. Files written: `CH4C.itp` (topology), `CH4C_cg.gro` (one molecule), `cg_groups.json` and `cg_index.json` (helper data). The zip already contains these, so re-running reproduces them.

### 3.2 What the mapping is

| Bead(s) | Type | Atoms (mol2 names) | Role |
|---|---|---|---|
| ETH, R1–R6, C1, C2 | N2a, SC4, SC3, TC2, C2 | steroid, from the Martini 3 cholesterol template | ETH replaces the cholesterol OH bead; covers C3 and the ether O |
| LK | C1 | C25 C27 C28 C29 | four-carbon linker |
| PY_A, PY_B, PY_C | TQ2 (+1), TC5, TC5 | pyridinium ring | charge on the bead holding N⁺ |
| P1_A, P1_B, P1_C | TC5 | first phenyl ring | three tiny beads per ring |
| VIN | SC4 | C12=C13 | alkene bead |
| CN | TN6a | C56≡N3 | nitrile |
| P2_A, P2_B, P2_C | TC5 | second phenyl ring | |
| PR_A, PR_B, PR_C | TC5, TC5, TN5a | pyridine ring | the acceptor bead holds the N |

Bead numbers 1–24 follow this order (ETH = 1, R1 = 2, …, C2 = 9, LK = 10, PY_A = 11, …, PR_C = 24). The bead types and bonded values are my starting choices. Rings are three tiny beads with 0.29 nm constraints (from my recollection of the Martini 3 benzene model; confirm it in the Martini 3 small-molecule paper). Linking bonds are at least 0.33 nm so bonded beads do not overlap. One improper dihedral at the vinyl bead keeps the E/Z configuration of your mol2.

### 3.3 Check the topology

```bash
python3 $S/check_cg_itp.py CH4C.itp $ROOT/ff/martini_v3.0.0.itp --charge 1
```
Expect: `24 beads`, `net charge +1.00`, `OK`. It confirms every bead type exists in your force-field file, every index is valid and nothing is disconnected.

### 3.4 Look at it

Open `CH4C_cg.gro` in VMD. Show the beads as spheres. You should see a compact steroid cluster, a short linker, and a flat rod-like head group with the nitrile bead sticking out from the middle.

### 3.5 Bead types to test later

| Choice | Starting value | Alternatives |
|---|---|---|
| Pyridinium N⁺ bead | TQ2 | TQ3 |
| Ether bead ETH | N2a | N3a |
| Nitrile bead CN | TN6a | TN5a |
| Pyridine N bead PR_C | TN5a | TN6a |
| Ring carbons | TC5 | TC4 for the pyridinium ring carbons |

To change a type, edit the `type` field in `$ROOT/params/chol4csn_cg_map.json` (`eth_type` for ETH) and repeat 3.1 and 3.3. You do not need to do this now; Part 9 uses it.

### 3.6 Topology header and the bromide structure

```bash
cd $ROOT/params
cat > topol_head.top <<'EOT'
#include "martini_v3.0.0.itp"
#include "martini_v3.0.0_ions_v1.itp"
#include "martini_v3.0.0_solvents_v1.itp"
#include "CH4C.itp"
EOT
cat > br.gro <<'EOT'
Br-
    1
    1ION     BR    1   0.000   0.000   0.000
   0.00000   0.00000   0.00000
EOT
```
The column positions in `br.gro` matter; copy it exactly. Atom names in a `.gro` file must match the topology or `grompp` stops: they do here (`BR` in residue `ION`; `W` in `W`; `S1`, `O2` in `DMS`).

---

## Part 4. Solvent boxes (build once, reuse)

### 4.1 Create the raw boxes

```bash
cd $ROOT/solvents
cp $ROOT/ff/martini_v3.0.0*.itp $ROOT/params/{CH4C.itp,topol_head.top} .
python3 $S/make_solvent_cg.py --n 10 --x 0.0 --out water_raw.gro
python3 $S/make_solvent_cg.py --n 10 --x 1.0 --out dmso_raw.gro
python3 $S/make_solvent_cg.py --n 10 --x 0.1 --out mix10_raw.gro
python3 $S/make_solvent_cg.py --n 10 --x 0.3 --out mix30_raw.gro
python3 $S/make_solvent_cg.py --n 10 --x 0.5 --out mix50_raw.gro
```
`x` is the DMSO mole fraction counting one W bead as 4 water molecules. Expect these counts: water 1000 W; DMSO 1000 molecules; mix10 692 W + 308 DMSO; mix30 368 W + 632 DMSO; mix50 200 W + 800 DMSO. Each box is 4.9 nm on a side.

### 4.2 Minimize and equilibrate each box

Run this block five times, with `NAME` set to `water`, `dmso`, `mix10`, `mix30`, `mix50`. Start with water and look at the result before doing the rest.

```bash
NAME=water
(cat topol_head.top; printf "\n[ system ]\n$NAME\n\n"; python3 $S/gro_to_molecules.py ${NAME}_raw.gro DMS=DMSO) > topol_$NAME.top
gmx grompp -f $MDP/em.mdp  -c ${NAME}_raw.gro -p topol_$NAME.top -o ${NAME}_em.tpr
gmx mdrun  -deffnm ${NAME}_em $MDRUN_OPTS
gmx grompp -f $MDP/eq1.mdp -c ${NAME}_em.gro  -p topol_$NAME.top -o ${NAME}_eq1.tpr
gmx mdrun  -deffnm ${NAME}_eq1 $MDRUN_OPTS
gmx grompp -f $MDP/eq2.mdp -c ${NAME}_eq1.gro -t ${NAME}_eq1.cpt -p topol_$NAME.top -o ${NAME}_eq2.tpr
gmx mdrun  -deffnm ${NAME}_eq2 $MDRUN_OPTS
cp ${NAME}_eq2.gro ${NAME}_box.gro
tail -n 1 ${NAME}_box.gro
```
What each step does: `em` removes overlaps (steepest-descent minimization); `eq1` is a 1 ns run with a 10 fs time step; `eq2` is a 5 ns run at 20 fs. Expect: `em` ends with "converged to Fmax < 1000" (or stops at the step limit with a small force); no LINCS warnings; the last line is the box size, a little different from 4.9 nm.

If `grompp` stops with a name mismatch or an unknown name, do not use `-maxwarn`; see Troubleshooting. (Residues named `DMS` map to the molecule `DMSO`; that is what `DMS=DMSO` does.)

### 4.3 Check the density

GROMACS' mass density is not meaningful here because Martini bead masses are nominal. Count molecules per volume instead: molecules ÷ (the three box edges multiplied). Expected number density:

- water: about 8.3 W beads per nm³ (33.3 waters per nm³ divided by 4)
- DMSO: about 8.4–8.5 molecules per nm³ (1.10 g/cm³)
- mixtures: between the two

With the 1000-molecule boxes this means box volume near 120 nm³ (edge 4.9–5.0 nm) for water and near 118 nm³ for DMSO. A box far from this means the topology or the run is wrong. Also open each mixture in VMD and confirm that water and DMSO are mixed, not separated into two slabs. The Martini 3 DMSO model in water/DMSO mixtures is less tested than in water, so state this when you report mixture results.

---

## Part 5. Build, equilibrate and validate any system

This part is a template. First use it for one molecule (5.5), later for dimers (Part 6) and the 24-molecule boxes (Part 7).

### 5.1 Choose the system

```bash
export SOLV=water                 # water | dmso | mix10 | mix30 | mix50
export NMOL=1 NION=1 L=6 SEED=1   # one molecule, one Br-, 6 nm box (for the first test)
export RUN=$ROOT/test/$SOLV
mkdir -p $RUN && cd $RUN
cp $ROOT/ff/martini_v3.0.0*.itp $ROOT/params/{CH4C.itp,CH4C_cg.gro,topol_head.top,br.gro,cg_groups.json} .
```
Later values: dimer `NMOL=2 NION=2 L=10`; aggregation `NMOL=24 NION=24 L=10`.

### 5.2 Place the molecules, bromide and solvent

```bash
gmx insert-molecules -ci CH4C_cg.gro -box $L $L $L -nmol $NMOL -scale 0.5 -try 2000 -rot xyz -seed $SEED -o solutes.gro
gmx insert-molecules -f solutes.gro -ci br.gro -nmol $NION -scale 0.5 -try 2000 -seed $SEED -o solutes_br.gro
gmx solvate -cp solutes_br.gro -cs $ROOT/solvents/${SOLV}_box.gro -radius 0.21 -o system.gro
```
Expect: `Added 1 molecules` (or `$NMOL`, `$NION`) for each `insert-molecules`, and `gmx solvate` reports the number of solvent molecules added. The number placed must equal what you asked for.

### 5.3 Topology and index file

```bash
(cat topol_head.top; printf "\n[ system ]\nChol4CSN in $SOLV\n\n"; python3 $S/gro_to_molecules.py system.gro ION=BR DMS=DMSO) > topol.top
tail -n 6 topol.top
python3 $S/make_ndx_cg.py cg_groups.json index.ndx --nmol $NMOL --nion $NION --gro system.gro
```
Expect: the `[ molecules ]` block lists `CH4C 1`, `BR 1`, then `W` and/or `DMSO` with their counts. The index file holds the groups used later (MOL1, MOL2, SOLUTE, STER, HEAD, NPLUS, NITRILE, RING1–4, ION, W, DMS, SOLUTE_BR).

### 5.4 Minimize and equilibrate

```bash
gmx grompp -f $MDP/em.mdp  -c system.gro -p topol.top -n index.ndx -o em.tpr
gmx mdrun  -deffnm em $MDRUN_OPTS
gmx grompp -f $MDP/eq1.mdp -c em.gro  -p topol.top -n index.ndx -o eq1.tpr
gmx mdrun  -deffnm eq1 $MDRUN_OPTS
gmx grompp -f $MDP/eq2.mdp -c eq1.gro -t eq1.cpt -p topol.top -n index.ndx -o eq2.tpr
gmx mdrun  -deffnm eq2 $MDRUN_OPTS
```
Expect: `grompp` reports the system is neutral (the +1 cations and −1 bromides cancel). If minimization stalls with very large forces, raise `-try` in 5.2 or lower `emstep` in `em.mdp` to 0.0005.

### 5.5 First validation: one molecule for 1 µs

Run this for `SOLV=water` first, then `SOLV=dmso`. At the end of 5.4 you are in `$ROOT/test/$SOLV`.

```bash
sed 's/nsteps .*/nsteps                   = 50000000/' $MDP/md_prod.mdp > md_test.mdp      # 1 us
gmx grompp -f md_test.mdp -c eq2.gro -p topol.top -n index.ndx -o md.tpr
gmx mdrun  -deffnm md $MDRUN_OPTS                      # on a CPU-only machine use a shorter run first, e.g. 100 ns
```
Then measure the shape and the E/Z configuration:
```bash
printf '[ ETH ]\n1\n[ PRC ]\n24\n' > ends.ndx
gmx distance -s md.tpr -f md.xtc -n ends.ndx -select 'com of group "ETH" plus com of group "PRC"' -oav ete.xvg
echo SOLUTE | gmx gyrate -s md.tpr -f md.xtc -n index.ndx -o rg.xvg
printf '[ ez ]\n16 18 17 19\n[ twist1 ]\n12 13 14 15\n[ twist2 ]\n20 21 22 23\n' > dih.ndx
gmx angle -f md.xtc -n dih.ndx -type dihedral -ov dih_time.xvg -od dih_dist.xvg
```
Plot `ete.xvg`, `rg.xvg` and `dih_time.xvg` (for example `python3 -c` with matplotlib, or `xmgrace` if installed).

**Pass if:**
- the end-to-end distance (ETH to PR_C) stays roughly constant at about 2.5–3 nm and the radius of gyration is stable (the molecule does not fold up);
- the `ez` dihedral stays within about 30° of its starting value (close to ±180°); a flip means the E/Z configuration was lost;
- the beads show no overlaps or LINCS errors in `md.log`.

If it fails, change bead types (3.5) or the bonded terms (9.1) and repeat. Do not start any PMF or aggregation run with an unstable molecule.

---

## Part 6. Dimer PMF (repeat for each solvent)

The coordinate is the distance between the centres of mass of the two cations. Box 10 nm, start at 4.5 nm, windows out to 4.4 nm. Do water first.

### 6.1 Starting dimers in three orientations

```bash
export SOLV=water NMOL=2 NION=2 L=10 SEED=1
for ORIENT in parallel antiparallel endtoend; do
  export RUN=$ROOT/dimer/$SOLV/$ORIENT; mkdir -p $RUN; cd $RUN
  cp $ROOT/ff/martini_v3.0.0*.itp $ROOT/params/{CH4C.itp,CH4C_cg.gro,topol_head.top,br.gro,cg_groups.json} .
  python3 $S/make_dimer_cg.py CH4C_cg.gro CH4C.itp solutes.gro --mode $ORIENT --sep 4.5 --box $L | tee dimer.log
done
```
Expect for each orientation: a line with the COM distance 4.500 nm and the closest bead contact, then `pull-group1-pbcatom = 12` and `pull-group2-pbcatom = 36`.

Then in each of the three folders do 5.2 (skip the first `insert-molecules`, since `solutes.gro` already exists: start from the bromide insertion), 5.3 and 5.4 with `NMOL=2 NION=2 L=10`:
```bash
cd $ROOT/dimer/water/parallel          # then antiparallel, then endtoend
export SOLV=water NMOL=2 NION=2 L=10 SEED=1
gmx insert-molecules -f solutes.gro -ci br.gro -nmol $NION -scale 0.5 -try 2000 -seed $SEED -o solutes_br.gro
gmx solvate -cp solutes_br.gro -cs $ROOT/solvents/${SOLV}_box.gro -radius 0.21 -o system.gro
(cat topol_head.top; printf "\n[ system ]\nChol4CSN dimer in $SOLV\n\n"; python3 $S/gro_to_molecules.py system.gro ION=BR DMS=DMSO) > topol.top
python3 $S/make_ndx_cg.py cg_groups.json index.ndx --nmol $NMOL --nion $NION --gro system.gro
gmx grompp -f $MDP/em.mdp  -c system.gro -p topol.top -n index.ndx -o em.tpr  && gmx mdrun -deffnm em  $MDRUN_OPTS
gmx grompp -f $MDP/eq1.mdp -c em.gro  -p topol.top -n index.ndx -o eq1.tpr && gmx mdrun -deffnm eq1 $MDRUN_OPTS
gmx grompp -f $MDP/eq2.mdp -c eq1.gro -t eq1.cpt -p topol.top -n index.ndx -o eq2.tpr && gmx mdrun -deffnm eq2 $MDRUN_OPTS
```

### 6.2 Check the start distance

```bash
gmx distance -s eq2.tpr -f eq2.gro -n index.ndx -select 'com of group "MOL1" plus com of group "MOL2"' -oav d0.xvg
tail -n 2 d0.xvg
```
Expect a value close to 4.5 nm. A big drift means the two molecules moved apart or together during equilibration; rebuild with a different `SEED`.

### 6.3 Steered pull (4.5 → 0.4 nm in 4.1 ns)

This generates the starting structures for all windows.
```bash
sed "s/PBC1/12/; s/PBC2/36/" $MDP/pull_steered.mdp > pull.mdp
gmx grompp -f pull.mdp -c eq2.gro -t eq2.cpt -p topol.top -n index.ndx -o pull.tpr
gmx mdrun  -deffnm pull -px pull_pullx.xvg -pf pull_pullf.xvg $MDRUN_OPTS
```
Plot the second column of `pull_pullx.xvg` (distance) against the first (time, ps). Expect it to fall from 4.5 nm and flatten where the molecules touch. That flat value is the closest approach. If it never flattens, the pull was too fast or the molecules slid past each other; this is normal for the end-to-end start, check the other two.

### 6.4 Create the umbrella windows

```bash
python3 $S/setup_windows_cg.py --pullx pull_pullx.xvg --xtc pull.xtc --tpr pull.tpr --ndx index.ndx --top topol.top \
        --pbc1 12 --pbc2 36 --ns 100 --outdir windows
head -n 8 windows/windows.csv
ls windows | head
```
Expect `43 windows written`. The schedule is 0.5–1.5 nm every 0.05 nm (spring constant 1000 kJ/mol/nm²), 1.6–3.0 nm every 0.1 nm (500), 3.2–4.4 nm every 0.2 nm (250). If the script lists windows with "no steered frame near the target", then the closest approach in 6.3 was above 0.5 nm; shift the first segment up, for example `--segments "0.65:1.5:0.05:1000,1.6:3.0:0.1:500,3.2:4.4:0.2:250"`.

`--ns 100` is the first pass (100 ns per window). Use `--ns 500` for the final result.

### 6.5 Run the windows

Serial, one after another (works anywhere):
```bash
cd $ROOT/dimer/water/parallel/windows
for w in win_*; do (cd $w && gmx mdrun -deffnm umb -px umb_pullx.xvg -pf umb_pullf.xvg $MDRUN_OPTS); done
```
On a cluster, run them as a job array. Save this as `umb.sh` and submit with `sbatch umb.sh` (edit the account and partition lines for your site):
```bash
#!/bin/bash
#SBATCH --job-name=umb
#SBATCH --array=0-42
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
source $HOME/chol4csn_m3/env.sh
cd $ROOT/dimer/water/parallel/windows/win_$(printf %03d $SLURM_ARRAY_TASK_ID)
gmx mdrun -deffnm umb -cpi umb.cpt -px umb_pullx.xvg -pf umb_pullf.xvg -ntmpi 1 -ntomp $SLURM_CPUS_PER_TASK -nb gpu
```
If a window is interrupted, run its `gmx mdrun` line again with `-cpi umb.cpt`; it resumes.

### 6.6 WHAM: join the windows into a PMF

```bash
cd $ROOT/dimer/water/parallel/windows
gmx wham -it tpr-files.dat -ix pullx-files.dat -o profile.xvg -hist histo.xvg -unit kJ -temp 300 \
         -b 20000 -bins 400 -nBootstrap 100 -bs-method b-hist -bsres bsres.xvg -bsprof bsprofs.xvg
```
`-b` is in ps: with 100 ns windows `-b 20000` discards the first 20 ns; with 500 ns windows use `-b 100000`. Plot `histo.xvg`: every histogram must overlap its neighbours. If there is a gap, add a window in the gap or raise the spring constant there, and rerun those windows.

### 6.7 Convert to a physical PMF and a binding free energy

```bash
python3 $S/pmf_analysis.py profile.xvg --bsprof bsprofs.xvg --plateau 3.8 4.4 --rb 1.0 1.5 2.0 --out pmf_physical.csv
python3 $S/plot_results.py pmf --histo histo.xvg --pmf pmf_physical.csv -o pmf.png
```
The script adds the geometric correction 2kT ln r, sets the zero of the PMF on the plateau (3.8–4.4 nm), integrates over the bound region up to each `--rb` (the cutoff between bound and unbound), and gives ΔG° = −kT ln(K·C°) with C° = 0.6022 nm⁻³ (the 1 M standard state). Report the well depth ± the bootstrap error as the headline number, and ΔG° with its dependence on `--rb`. With the reaction-field electrostatics cut off at 1.1 nm, the PMF should be flat well before 3.8 nm; if it still slopes there, extend the windows.

### 6.8 Convergence

Use the same window length you used in 6.6. For 500 ns windows:
```bash
gmx wham -it tpr-files.dat -ix pullx-files.dat -o profile_h1.xvg -unit kJ -temp 300 -b 100000 -e 300000 -bins 400
gmx wham -it tpr-files.dat -ix pullx-files.dat -o profile_h2.xvg -unit kJ -temp 300 -b 300000 -e 500000 -bins 400
python3 $S/pmf_analysis.py profile_h1.xvg --plateau 3.8 4.4 --out pmf_h1.csv
python3 $S/pmf_analysis.py profile_h2.xvg --plateau 3.8 4.4 --out pmf_h2.csv
python3 $S/plot_results.py pmf --histo histo.xvg --pmf pmf_h1.csv --pmf2 pmf_h2.csv -o pmf_halves.png
```
(With 100 ns windows use `-b 20000 -e 60000` and `-b 60000 -e 100000`.) The two halves must agree within the bootstrap error (well depths within about 2 kJ/mol). If not, run longer windows (`gmx convert-tpr -s umb.tpr -extend <ps> -o umb.tpr`, then `mdrun -cpi`). Also compare the three orientations: a spread above about 5 kJ/mol means the distance-only PMF is limited by how the molecules are oriented relative to each other, and you should say so when reporting.

### 6.9 Other solvents

Repeat 6.1–6.8 with `SOLV=dmso`, then `mix10`, `mix30`, `mix50`. Compare well depth, minimum position, barrier height and ΔG° across solvents.

---

## Part 7. Aggregation simulations

### 7.1 Build three independent boxes per solvent (water first, then DMSO)

For each `SEED` in 1, 2, 3, run 5.1–5.4 with these values:
```bash
export SOLV=water NMOL=24 NION=24 L=10
export SEED=1
export RUN=$ROOT/agg/$SOLV/rep$SEED
mkdir -p $RUN && cd $RUN
cp $ROOT/ff/martini_v3.0.0*.itp $ROOT/params/{CH4C.itp,CH4C_cg.gro,topol_head.top,br.gro,cg_groups.json} .
# then 5.2 (both insert-molecules and solvate), 5.3 and 5.4
```
Concentration: c = NMOL / (0.6022 × L³) mol/L, so 24 molecules in a 10 nm box is about 40 mM. If you know the experimental concentration, choose NMOL and L to match it. Expect about 9,000 beads in water and 17,000 in DMSO.

### 7.2 Production run, 5 µs each

```bash
cd $ROOT/agg/water/rep1
gmx grompp -f $MDP/md_prod.mdp -c eq2.gro -t eq2.cpt -p topol.top -n index.ndx -o md.tpr
gmx mdrun  -deffnm md -cpt 15 $MDRUN_OPTS
# to resume after an interruption:  gmx mdrun -deffnm md -cpi md.cpt $MDRUN_OPTS
```
Benchmark first (about a minute) so you know how long it will take:
```bash
gmx mdrun -deffnm bench -s md.tpr -nsteps 50000 -resethway -noconfout $MDRUN_OPTS
grep Performance bench.log           # first number is ns/day; divide by 1000 for µs/day
```

### 7.3 Extend if the aggregation is still developing

Aggregation is usually slower in DMSO.
```bash
gmx convert-tpr -s md.tpr -extend 2500000 -o md.tpr        # +2.5 us (the value is in ps)
gmx mdrun -deffnm md -cpi md.cpt $MDRUN_OPTS
```

### 7.4 Solvent-swap test (2 µs each way)

This tests whether an aggregate formed in one solvent survives in the other.
```bash
cd $ROOT/agg/water/rep1
echo SOLUTE_BR | gmx trjconv -s md.tpr -f md.xtc -n index.ndx -pbc mol -dump 5000000 -o last_solute.gro      # time in ps
gmx editconf -f last_solute.gro -box 10 10 10 -o solutes_br.gro
export SOLV=dmso NMOL=24 NION=24 L=10 SEED=1 RUN=$ROOT/swap/water_to_dmso
mkdir -p $RUN && cd $RUN
cp $ROOT/agg/water/rep1/solutes_br.gro .
cp $ROOT/ff/martini_v3.0.0*.itp $ROOT/params/{CH4C.itp,CH4C_cg.gro,topol_head.top,br.gro,cg_groups.json} .
# then: solvate (5.2 last line, with the DMSO box), 5.3, 5.4, and a 2 us run:
sed 's/nsteps .*/nsteps                   = 100000000/' $MDP/md_prod.mdp > md_swap.mdp
```
Then use `md_swap.mdp` in place of `md_prod.mdp` as in 7.2. Do the reverse (DMSO to water) too. If the cluster is cut across the box edge in `last_solute.gro`, first apply `gmx trjconv -pbc cluster`.

---

## Part 8. Analysis

### 8.1 Aggregation metrics

```bash
cd $ROOT/agg/water/rep1
echo SOLUTE_BR | gmx trjconv -s md.tpr -f md.xtc -n index.ndx -pbc mol -o solute_whole.xtc
echo SOLUTE_BR | gmx trjconv -s md.tpr -f md.xtc -n index.ndx -pbc mol -dump 0 -o solute.gro
python3 $S/analyze_agg_cg.py --top solute.gro --traj solute_whole.xtc --groups cg_groups.json --itp CH4C.itp \
        --nmol 24 --nion 24 --cutoff 0.6 --stride 5 --out agg_rep1.csv
python3 $S/summarize_reps.py $ROOT/agg/water/rep{1,2,3}/agg_rep*.csv --last-ns 500 --nmol 24
python3 $S/plot_results.py agg $ROOT/agg/water/rep{1,2,3}/agg_rep*.csv -o agg_water.png
```
Reading `.xtc` needs `pip install MDAnalysis`. To avoid installing it, change the output name of the first `trjconv` line to `solute_whole.gro` and pass that as `--traj`; the script reads multi-frame `.gro` files itself. Run the first three commands in each replicate folder (name the csv `agg_rep1.csv`, `agg_rep2.csv`, `agg_rep3.csv`), then the last two once.

Columns in the csv: time, number of clusters, monomers, largest cluster, weight-average cluster size, radius of gyration and shape factor κ² of the largest cluster, stacked ring pairs, Br⁻ per cation. A cluster is a set of molecules joined by any bead–bead distance below `--cutoff` (default 0.6 nm). Repeat with `--cutoff 0.5` and `0.7`; the conclusions should not change. `--last-ns` is in ns of simulation time.

### 8.2 Radial distribution functions and coordination numbers

Use the second half of the run (for 5 µs, `b=2500000` ps):
```bash
b=2500000
gmx rdf -s md.tpr -f md.xtc -n index.ndx -b $b -ref NPLUS   -sel ION  -bin 0.01 -o rdf_Nplus_Br.xvg
gmx rdf -s md.tpr -f md.xtc -n index.ndx -b $b -ref HEAD    -sel HEAD -bin 0.01 -o rdf_head_head.xvg
gmx rdf -s md.tpr -f md.xtc -n index.ndx -b $b -ref STER    -sel STER -bin 0.01 -o rdf_ster_ster.xvg
gmx rdf -s md.tpr -f md.xtc -n index.ndx -b $b -ref NPLUS   -sel W    -bin 0.01 -cn -o rdf_Nplus_W.xvg       # water systems
gmx rdf -s md.tpr -f md.xtc -n index.ndx -b $b -ref NPLUS   -sel DMS  -bin 0.01 -cn -o rdf_Nplus_DMS.xvg     # DMSO systems
gmx rdf -s md.tpr -f md.xtc -n index.ndx -b $b -ref NITRILE -sel W    -bin 0.01 -cn -o rdf_CN_W.xvg
```
Martini has no hydrogen bonds, so use these coordination numbers (the `-cn` output) instead of hydrogen-bond counts. Run the `W` lines only in systems that contain water, and the `DMS` line only in systems that contain DMSO.

### 8.3 Buried surface (optional)

```bash
echo "??? ??? 0.264" > vdwradii.dat
echo SOLUTE | gmx sasa -s md.tpr -f md.xtc -n index.ndx -probe 0.191 -ndots 24 -o sasa.xvg
```
Martini beads need an explicit radius file. Test on one molecule first because the radius convention changes the absolute value; report the trend against cluster size, not the absolute number.

### 8.4 What to report

PMF well depth ± error, minimum position, barrier and ΔG° per solvent with the orientation spread; aggregation mean ± standard error over three replicates (cluster size, time to first dimer, ion pairing, stacking); whether the PMF ordering across solvents matches the aggregation ordering; the bead-type choices; box sizes; cluster cutoff; plateau window; and the raw simulation time (do not rescale it to "effective" time).

---

## Part 9. Controls and later improvements

### 9.1 Refine the bonded terms from an all-atom run (optional, recommended later)

Your first-pass bonded values come from the mapped geometry. To base them on real fluctuations, run one molecule plus one Br⁻ in water for 50–100 ns with the all-atom CHARMM36/CGenFF protocol (the earlier document `Chol4CSN_GROMACS_runbook.md`, Parts A, B and D), then:
```bash
cd $ROOT/aa_ref
# choose the group containing only the 130 atoms of the molecule
gmx trjconv -s md.tpr -f md.xtc -n index.ndx -pbc mol -o aa_solute.gro
cd $ROOT/params
python3 $S/cg_bonded_from_aa.py --mol2 Chol4CSN_clean.mol2 --map chol4csn_cg_map.json --index cg_index.json \
        --gro $ROOT/aa_ref/aa_solute.gro --out bonded_overrides.json
python3 $S/cg_build.py Chol4CSN_clean.mol2 chol4csn_cg_map.json --outdir . --overrides bonded_overrides.json
python3 $S/check_cg_itp.py CH4C.itp $ROOT/ff/martini_v3.0.0.itp --charge 1
```
The script prints each bond and angle with its mean and spread, proposes force constants (k = kT/σ²), keeps the Martini bond lengths, and writes a torsion for each ring–ring bond whose distribution has one peak. Then repeat Part 5.5 with the new `CH4C.itp`. You can also compare the single-molecule end-to-end distance and radius of gyration with the all-atom run (atoms 28 29 and 23 24 in the mol2 numbering correspond to the ETH and PR_C ends):
```bash
gmx distance -s aa.tpr -f aa.xtc -select 'com of atomnr 28 29 plus com of atomnr 23 24' -oav ete_aa.xvg
```
Agreement within about 10–15% is a pass.

### 9.2 Bead-type sensitivity (do this before you write anything up)

Repeat the water dimer PMF (parallel orientation, 100 ns windows) with one alternative from 3.5, for example TQ3 for the N⁺ bead. Report both curves. A change in well depth larger than the run-to-run error is a modelling uncertainty you must quote.

### 9.3 Calibration against data

If you have an experimental critical aggregation concentration, aggregate size, or any all-atom dimer free energy in water, compare it with the CG result. A CG well that is more than about 5 kJ/mol deeper or shallower than the reference is a reason to adjust bead types (3.5).

### 9.4 Other checks

- Cluster definition: three cutoffs (8.1). Random seeds: three replicates (7.1).
- Look at the trend across water, mixtures and DMSO, not at a single absolute ΔG; trends are more robust than absolute values.

---

## Part 10. Cost, order and limits

**Cost.** Windows = 43 × (window length in µs) ÷ P days per orientation per solvent; aggregation = 3 × 5 µs ÷ P days per solvent, where P is your benchmark in µs/day from 7.2. At P = 5 µs/day, 500 ns windows cost about 4 GPU-days per orientation and aggregation about 3 GPU-days per solvent. I have not benchmarked this system, so P is a placeholder; use your own number.

**Order of work.**
1. Parts 1–5 completely, water and DMSO single-molecule runs included.
2. Water: first-pass PMF (100 ns windows, one orientation), then the bead-type sensitivity (9.2), then the full PMF, then water aggregation.
3. DMSO: PMF, then aggregation.
4. Mixtures, as time allows.
5. Analysis and comparison (Part 8).

**Limits you must state when reporting.**
- Martini 3 uses one dielectric constant (ε_r = 15 with reaction-field) in every solvent. The real change from water (ε ≈ 78) to DMSO (ε ≈ 47) is not modelled, so Br⁻ pairing with the cation in DMSO is probably under-represented. Solvent differences come only from the bead types.
- The Martini 3 DMSO model and the Martini 3 bromide bead were not parametrized together with your molecule. Check water/DMSO miscibility before using the mixtures.
- Your molecule's bead types and bonded terms are estimates until you have tested them (5.5, 9.1, 9.2).
- The PMF is a coarse-grained free energy. Use it to compare solvents, not as an absolute binding free energy.
- Martini dynamics are faster than real ones; do not convert simulation time into real time.

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| `gmx: command not found` | GROMACS not loaded in this terminal | `source ~/chol4csn_m3/env.sh` |
| `grompp`: atom name does not match | `.gro` and topology names differ | use the files exactly as in 3.6, 4.1, 5.2; do not use `-maxwarn` |
| `grompp`: unknown atomtype | a bead type in `CH4C.itp` is not in the force-field file | `check_cg_itp.py` names it |
| `grompp`: total charge is not zero | `NION` differs from `NMOL` | set `NION = NMOL` and rebuild |
| `grompp`: cannot find include file | missing `.itp` in the run folder | repeat the `cp` line of 5.1 |
| LINCS warnings | overlaps or a large first step | rerun `em.mdp`, then `eq1.mdp` (10 fs) before `eq2.mdp` |
| `mdrun`: "no GPU"/option `-nb gpu` unknown | CPU-only GROMACS | remove `-nb gpu` in `env.sh` |
| steroid distorts | edited cholesterol part | regenerate `CH4C.itp` with 3.1 |
| pull distance exceeds half the box | box too small | use a box more than twice the largest pull distance |
| gaps between umbrella histograms | windows too sparse | add windows or raise the spring constant |
| clusters broken across the box | trajectory not made whole | use `-pbc mol` in `trjconv` |
| solvent box demixes | untested water/DMSO mixture | equilibrate longer and report it |
| `insert-molecules` places fewer molecules | box crowded | raise `-try`, or use a bigger box |

---

## Appendix: the five mdp files

Save each block under the name shown in `$MDP`; they are also in the zip. They are Martini 3 settings: reaction-field electrostatics (ε_r = 15), cutoff 1.1 nm, V-rescale thermostat at 300 K, 20 fs time step in production.

### em.mdp
```
; Martini 3 energy minimisation. FLEXIBLE turns the ring and steroid constraints into stiff bonds.
define                   = -DFLEXIBLE
integrator               = steep
emtol                    = 1000
emstep                   = 0.001
nsteps                   = 20000
cutoff-scheme            = Verlet
verlet-buffer-tolerance  = 0.005
coulombtype              = reaction-field
rcoulomb                 = 1.1
epsilon-r                = 15
epsilon-rf               = 0
vdwtype                  = cutoff
vdw-modifier             = potential-shift-verlet
rvdw                     = 1.1
constraints              = none
constraint-algorithm     = LINCS
lincs-order              = 12
lincs-iter               = 2
pbc                      = xyz
```

### eq1.mdp
```
; Martini 3 equilibration 1: 1 ns at dt = 10 fs, C-rescale barostat (also used for the pure-solvent boxes)
integrator               = md
dt                       = 0.010
nsteps                   = 100000
nstlog                   = 5000
nstenergy                = 1000
nstxout-compressed       = 10000
gen-vel                  = yes
gen-temp                 = 300
gen-seed                 = -1
continuation             = no
cutoff-scheme            = Verlet
verlet-buffer-tolerance  = 0.005
coulombtype              = reaction-field
rcoulomb                 = 1.1
epsilon-r                = 15
epsilon-rf               = 0
vdwtype                  = cutoff
vdw-modifier             = potential-shift-verlet
rvdw                     = 1.1
constraints              = none
constraint-algorithm     = LINCS
lincs-order              = 12
lincs-iter               = 2
pbc                      = xyz
tcoupl                   = v-rescale
tc-grps                  = System
tau-t                    = 1.0
ref-t                    = 300
pcoupl                   = c-rescale
pcoupltype               = isotropic
tau-p                    = 4.0
ref-p                    = 1.0
compressibility          = 3e-4
```

### eq2.mdp
```
; Martini 3 equilibration 2: 5 ns at dt = 20 fs
integrator               = md
dt                       = 0.020
nsteps                   = 250000
nstlog                   = 5000
nstenergy                = 1000
nstxout-compressed       = 5000
continuation             = yes
cutoff-scheme            = Verlet
verlet-buffer-tolerance  = 0.005
coulombtype              = reaction-field
rcoulomb                 = 1.1
epsilon-r                = 15
epsilon-rf               = 0
vdwtype                  = cutoff
vdw-modifier             = potential-shift-verlet
rvdw                     = 1.1
constraints              = none
constraint-algorithm     = LINCS
lincs-order              = 12
lincs-iter               = 2
pbc                      = xyz
tcoupl                   = v-rescale
tc-grps                  = System
tau-t                    = 1.0
ref-t                    = 300
pcoupl                   = c-rescale
pcoupltype               = isotropic
tau-p                    = 4.0
ref-p                    = 1.0
compressibility          = 3e-4
```

### md_prod.mdp
```
; Martini 3 production, 24-molecule box.  nsteps = 250000000 is 5 us at dt = 20 fs.
integrator               = md
dt                       = 0.020
nsteps                   = 250000000
nstlog                   = 50000
nstenergy                = 5000
nstxout-compressed       = 5000             ; every 100 ps, whole system
continuation             = no
gen-vel                  = yes
gen-temp                 = 300
gen-seed                 = -1
cutoff-scheme            = Verlet
verlet-buffer-tolerance  = 0.005
coulombtype              = reaction-field
rcoulomb                 = 1.1
epsilon-r                = 15
epsilon-rf               = 0
vdwtype                  = cutoff
vdw-modifier             = potential-shift-verlet
rvdw                     = 1.1
constraints              = none
constraint-algorithm     = LINCS
lincs-order              = 12
lincs-iter               = 2
pbc                      = xyz
tcoupl                   = v-rescale
tc-grps                  = System
tau-t                    = 1.0
ref-t                    = 300
pcoupl                   = parrinello-rahman
pcoupltype               = isotropic
tau-p                    = 12.0
ref-p                    = 1.0
compressibility          = 3e-4
```

### pull_steered.mdp
```
; Martini 3 steered pull: the two solutes are drawn from 4.5 nm to 0.4 nm COM distance in 4.1 ns
integrator               = md
dt                       = 0.020
nsteps                   = 205000
nstlog                   = 5000
nstenergy                = 1000
nstxout-compressed       = 100             ; every 2 ps, whole system (windows start from these frames)
continuation             = yes
cutoff-scheme            = Verlet
verlet-buffer-tolerance  = 0.005
coulombtype              = reaction-field
rcoulomb                 = 1.1
epsilon-r                = 15
epsilon-rf               = 0
vdwtype                  = cutoff
vdw-modifier             = potential-shift-verlet
rvdw                     = 1.1
constraints              = none
constraint-algorithm     = LINCS
lincs-order              = 12
lincs-iter               = 2
pbc                      = xyz
tcoupl                   = v-rescale
tc-grps                  = System
tau-t                    = 1.0
ref-t                    = 300
pcoupl                   = parrinello-rahman
pcoupltype               = isotropic
tau-p                    = 12.0
ref-p                    = 1.0
compressibility          = 3e-4
pull                     = yes
pull-ngroups             = 2
pull-ncoords             = 1
pull-group1-name         = MOL1
pull-group2-name         = MOL2
pull-group1-pbcatom      = PBC1             ; replace with the numbers printed by make_dimer_cg.py
pull-group2-pbcatom      = PBC2
pull-coord1-type         = umbrella
pull-coord1-geometry     = distance
pull-coord1-groups       = 1 2
pull-coord1-dim          = Y Y Y
pull-coord1-start        = yes
pull-coord1-rate         = -0.001           ; nm/ps, negative = the groups approach each other
pull-coord1-k            = 1000
pull-nstxout             = 50
pull-nstfout             = 50
```
