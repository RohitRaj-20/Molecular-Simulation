# Third-party material

- `scripts/cg_build.py` embeds the Martini 3 cholesterol (CHOL) parameters from `martini_v3.0.0_sterols_v1.itp` in the Martini Force Field Initiative repository `M3-Lipid-Parameters`, which is licensed under the Apache License 2.0. Those parameters are the work of the Martini 3 cholesterol authors (Borges-Araújo et al., 2023); please cite that paper if you use this model.
- The Martini 3 force-field files (`martini_v3.0.0.itp`, ions, solvents) are not included. They are downloaded by `ff/get_martini_ff.sh` from the same repository and keep their own licence.
- Please also cite the Martini 3 force-field paper (Souza et al., Nature Methods, 2021) and GROMACS.
