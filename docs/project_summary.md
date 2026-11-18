# Project summary

**Title:** Self-Aggregation of a Cholesterol–α-Cyanostilbene Cation in Water and DMSO Using Martini Coarse-Grained Simulations

**Description:** This project studies how a cholesterol-linked, positively charged α-cyanostilbene molecule (with bromide as the counter-ion) self-aggregates. Martini coarse-grained simulations of the whole molecule are run in GROMACS, first in water and then in DMSO, with water/DMSO mixtures also included. Cluster sizes, ring stacking and bromide positions are analysed, and the PMF of two molecules approaching each other is calculated in each solvent to explain the differences in aggregation.

**Challenge:** No Martini parameters exist for this molecule, so the model had to be built by mapping its 130 atoms onto 24 beads. Aggregation is slow and random, so several independent runs are needed per solvent. Martini also uses the same dielectric constant for every solvent, so water and DMSO differ only through their bead types. Solvent trends are compared instead of single exact values, and an alternative set of bead types is tested to check whether the results change.
