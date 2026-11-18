#!/usr/bin/env bash
# Rebuilds the Martini model of Chol4CSN and checks that it matches the stored files.
# Needs only python3 with numpy and scipy. If ff/ holds the Martini files, the topology is also checked.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; S=$ROOT/scripts; T=$(mktemp -d)
python3 "$S/cg_build.py" "$ROOT/params/Chol4CSN_clean.mol2" "$ROOT/params/chol4csn_cg_map.json" --outdir "$T" > "$T/build.log"
head -n 2 "$T/build.log"
diff -q "$T/CH4C.itp" "$ROOT/params/CH4C.itp" && echo "CH4C.itp reproduces"
diff -q "$T/CH4C_cg.gro" "$ROOT/params/CH4C_cg.gro" && echo "CH4C_cg.gro reproduces"
if [ -f "$ROOT/ff/martini_v3.0.0.itp" ]; then
  python3 "$S/check_cg_itp.py" "$ROOT/params/CH4C.itp" "$ROOT/ff/martini_v3.0.0.itp" --charge 1 > "$T/check.log"
  tail -n 1 "$T/check.log"
else
  echo "skipped topology check: run  bash ff/get_martini_ff.sh  first"
fi
python3 "$S/make_solvent_cg.py" --n 10 --x 0.3 --out "$T/mix.gro" > "$T/solv.log"
head -n 1 "$T/solv.log"
python3 "$S/make_dimer_cg.py" "$ROOT/params/CH4C_cg.gro" "$ROOT/params/CH4C.itp" "$T/dimer.gro" --mode parallel > "$T/dim.log"
head -n 1 "$T/dim.log"
rm -rf "$T"; echo "smoke test passed"
