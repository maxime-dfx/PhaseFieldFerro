#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=PhaseFieldFerro
#SBATCH --output=%x-%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00

set -euo pipefail

# Définition des variables par défaut
VARIANT="${1:-pardiso}" # pardiso, release, tracy, ou pardiso_tracy
CONFIG="${2:-$HOME/Min/input/config.toml}"

# Gestion de l'historique d'attachement des threads
if [[ "$VARIANT" == *"tracy"* ]]; then
    export OMP_PROC_BIND=close
else
    export OMP_PROC_BIND=spread
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# On remonte d'un dossier car common_env.sh est maintenant dans setup/
source "$SCRIPT_DIR/../setup/common_env.sh" "$VARIANT"

echo "Démarrage de la simulation ($VARIANT) sur $OMP_NUM_THREADS coeurs..."
cd "$MIN_BUILD_DIR"
"$MIN_BINARY" "$CONFIG"
echo "Simulation terminée."