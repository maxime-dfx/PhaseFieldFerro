#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=PhaseFieldFerro_pardiso_tracy
#SBATCH --output=%x-%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/common_env.sh" pardiso_tracy

CONFIG="${1:-$HOME/Min/input/config.toml}"

echo "Job SLURM $SLURM_JOB_ID sur $(hostname), $OMP_NUM_THREADS coeurs alloués (PARDISO + Tracy)."
cd "$MIN_BUILD_DIR"
"$MIN_BINARY" "$CONFIG"
