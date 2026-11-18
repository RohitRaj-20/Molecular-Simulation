# Source this in every new terminal:  source ~/chol4csn_m3/env.sh
export ROOT=$HOME/chol4csn_m3
export S=$ROOT/scripts
export MDP=$ROOT/mdp
[ -f $ROOT/venv/bin/activate ] && source $ROOT/venv/bin/activate
[ -f $HOME/gromacs/bin/GMXRC ] && source $HOME/gromacs/bin/GMXRC
# mdrun options: set -ntomp to your CPU cores (type: nproc). Remove "-nb gpu" if you have no GPU.
export MDRUN_OPTS="-ntmpi 1 -ntomp 8 -nb gpu"
