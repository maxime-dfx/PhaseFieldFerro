#!/usr/bin/env bash
# =============================================================================
# configure.sh — Configure CMake pour PhaseFieldFerro
#
# Ce script permet de générer rapidement des dossiers de build dédiés
# selon les combinaisons de solveurs et de profilage choisies.
#
# Options disponibles :
#   --pardiso  : Active Intel MKL PARDISO (solveur direct)
#   --petsc    : Active PETSc & Hypre (solveur itératif AMG/ILU/GMRES)
#   --tracy    : Active le profilage Tracy
#
# Exemples d'usage :
#   ./configure.sh                           -> build_standard (sans options)
#   ./configure.sh --pardiso --petsc         -> build_pardiso_petsc (PROD max perf)
#   ./configure.sh --tracy                   -> build_tracy
#
# Pour surcharger un chemin de dépendance (ex: Gmsh) :
#   ./configure.sh --pardiso -- -DGMSH_ROOT=/autre/chemin
# =============================================================================
set -euo pipefail

# Détection automatique du dossier racine (là où se trouve ce script)
SRC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

# Valeurs par défaut
USE_PARDISO="OFF"
USE_PETSC="OFF"
ENABLE_PROFILING="OFF"
EXTRA_ARGS=()

# Parsing des arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --pardiso) USE_PARDISO="ON"; shift ;;
        --petsc)   USE_PETSC="ON"; shift ;;
        --tracy)   ENABLE_PROFILING="ON"; shift ;;
        --)        shift; EXTRA_ARGS+=("$@"); break ;;
        *)         echo "Argument inconnu: $1" >&2; exit 1 ;;
    esac
done

# Construction dynamique du nom du dossier de build
DIR_SUFFIX=""
[[ "$USE_PARDISO" == "ON" ]] && DIR_SUFFIX="${DIR_SUFFIX}_pardiso"
[[ "$USE_PETSC" == "ON" ]]   && DIR_SUFFIX="${DIR_SUFFIX}_petsc"
[[ "$ENABLE_PROFILING" == "ON" ]] && DIR_SUFFIX="${DIR_SUFFIX}_tracy"

# Détermination du répertoire cible
if [[ -z "$DIR_SUFFIX" ]]; then
    BUILD_DIR="$HOME/Min/build_standard"
else
    BUILD_DIR="$HOME/Min/build${DIR_SUFFIX}"
fi

echo "=== Configuration CMake ==="
echo "Dossier source : $SRC_DIR"
echo "Dossier cible  : $BUILD_DIR"
echo "PARDISO        : $USE_PARDISO"
echo "PETSc          : $USE_PETSC"
echo "Profilage Tracy: $ENABLE_PROFILING"
[[ ${#EXTRA_ARGS[@]} -gt 0 ]] && echo "Args suppl.    : ${EXTRA_ARGS[*]}"
echo "==========================="

# Création et déplacement dans le dossier cible
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

# Lancement de CMake
cmake "$SRC_DIR" \
    -DCMAKE_BUILD_TYPE=Release \
    -DUSE_PARDISO="$USE_PARDISO" \
    -DUSE_PETSC="$USE_PETSC" \
    -DENABLE_PROFILING="$ENABLE_PROFILING" \
    "${EXTRA_ARGS[@]}"

echo
echo "Configuration terminée. Pour compiler :"
echo "  cd $BUILD_DIR && make -j\$(nproc)"