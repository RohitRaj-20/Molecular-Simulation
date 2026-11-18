#!/usr/bin/env bash
# Download the three Martini 3 force-field files this project needs into ff/.
set -euo pipefail
cd "$(dirname "$0")"
rm -rf M3-Lipid-Parameters
git clone --depth 1 https://github.com/Martini-Force-Field-Initiative/M3-Lipid-Parameters.git
cp M3-Lipid-Parameters/ITPs/martini_v3.0.0.itp \
   M3-Lipid-Parameters/ITPs/martini_v3.0.0_ions_v1.itp \
   M3-Lipid-Parameters/ITPs/martini_v3.0.0_solvents_v1.itp .
rm -rf M3-Lipid-Parameters
ls -1 *.itp
grep -n -A4 "BR" martini_v3.0.0_ions_v1.itp | head -n 8
grep -n -A6 "DMSO" martini_v3.0.0_solvents_v1.itp | head -n 10
