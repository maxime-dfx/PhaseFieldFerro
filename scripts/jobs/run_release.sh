#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=PhaseFieldFerro-release
#SBATCH --output=%x-%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#
# Build standard : SimplicialLDLT/ConjugateGradient (Eigen pur, sans MKL).

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export OMP_PROC_BIND=spread   # historique : spread ici, close dans les variantes tracy
source "$SCRIPT_DIR/common_env.sh" release

CONFIG="${1:-$HOME/Min/input/config.toml}"

echo "Démarrage de la simulation sur $OMP_NUM_THREADS coeurs (build standard, sans PARDISO)..."
cd "$MIN_BUILD_DIR"
"$MIN_BINARY" "$CONFIG"
echo "Simulation terminée."
