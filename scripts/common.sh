# Shared settings. Override MOCKCHECK_SCRATCH for your cluster.
export MOCKCHECK_SCRATCH="${MOCKCHECK_SCRATCH:-/scratch/$USER/mockcheck}"
export MC_ENV="$MOCKCHECK_SCRATCH/env"
export MC_DB="$MOCKCHECK_SCRATCH/db"
export MC_READS="$MOCKCHECK_SCRATCH/reads"
export MC_OUT="$MOCKCHECK_SCRATCH/out"
mc_activate() { module load miniforge3/26.7.2-0 2>/dev/null; eval "$(conda shell.bash hook)"; conda activate "$MC_ENV"; }
