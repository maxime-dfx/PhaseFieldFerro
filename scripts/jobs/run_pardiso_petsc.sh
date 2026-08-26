#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=PhaseFieldFerro_pardiso_petsc
#SBATCH --output=%x-%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00

set -eo pipefail

SCRIPT_DIR="$HOME/Min/scripts/setup" 
source "$SCRIPT_DIR/common_env.sh" pardiso_petsc

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# --- SOLUTIONS ANTI-CRASH INTERNE MPI ---
export FI_PROVIDER=tcp
export I_MPI_FABRICS=shm
# ----------------------------------------

CONFIG="${1:-$HOME/Min/input/hysteresis.toml}"

echo "Job SLURM $SLURM_JOB_ID sur $(hostname), $OMP_NUM_THREADS coeurs alloués."
cd "$MIN_BUILD_DIR"

# Lancement DIRECT
"$MIN_BINARY" "$CONFIG"