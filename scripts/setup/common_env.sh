#!/usr/bin/env bash
# =============================================================================
# common_env.sh — setup d'environnement partagé par tous les scripts SLURM
# (run.sh, run_release.sh, run_pardiso_tracy.sh, run_tracy_live.sh,
# run_bench.sh, bench_scaling.sh).
#
# But : un seul endroit pour le setup MKL/OpenMP/modules, au lieu de le
# dupliquer (et de le désynchroniser) dans chaque script SLURM.
#
# Usage (à sourcer, pas à exécuter) :
#   source "$(dirname "${BASH_SOURCE[0]}")/common_env.sh" <build_variant>
#
# <build_variant> ∈ {release, pardiso, tracy, pardiso_tracy}
#   release        -> build_release        (Eigen pur, sans MKL)
#   pardiso        -> build                (PARDISO/MKL, prod)
#   tracy          -> build_tracy          (Eigen pur + Tracy)
#   pardiso_tracy  -> build_pardiso_tracy  (PARDISO/MKL + Tracy)
#
# Après le source, les variables suivantes sont disponibles :
#   $MIN_BUILD_DIR   chemin absolu du dossier de build correspondant
#   $MIN_BINARY      chemin absolu de l'exécutable PhaseFieldFerro
#
# Variables SLURM/OMP déjà exportées : OMP_NUM_THREADS (= SLURM_CPUS_PER_TASK
# si défini), OMP_PROC_BIND=close, OMP_PLACES=cores. bench_scaling.sh peut
# les écraser ensuite pour son propre balayage — c'est voulu, dernier mot à
# l'appelant.
# =============================================================================
set -euo pipefail

MIN_VARIANT="${1:?Usage: source common_env.sh <release|pardiso|tracy|pardiso_tracy>}"

source /etc/profile.d/z99-modules.sh
module purge

case "$MIN_VARIANT" in
    release)
        MIN_BUILD_DIR="$HOME/Min/build_release"
        MIN_USE_MKL=0
        ;;
    pardiso)
        MIN_BUILD_DIR="$HOME/Min/build"
        MIN_USE_MKL=1
        ;;
    tracy)
        MIN_BUILD_DIR="$HOME/Min/build_tracy"
        MIN_USE_MKL=0
        ;;
    pardiso_tracy)
        MIN_BUILD_DIR="$HOME/Min/build_pardiso_tracy"
        MIN_USE_MKL=1
        ;;
    pardiso_petsc) 
        MIN_BUILD_DIR="$HOME/Min/build_pardiso_petsc"
        MIN_USE_MKL=1
        ;;
    *)
        echo "[common_env] variante inconnue: $MIN_VARIANT (attendu: release|pardiso|tracy|pardiso_tracy|pardiso_petsc)" >&2
        return 1 2>/dev/null || exit 1
        ;;
esac

if [[ "$MIN_USE_MKL" -eq 1 ]]; then
    module load mkl/latest
    # MKL : un seul runtime OpenMP dans le process (evite le conflit
    # libgomp / libiomp5), threading gere par les #pragma omp du code.
    export MKL_THREADING_LAYER=GNU
    export MKL_INTERFACE_LAYER=LP64
    export MKL_DYNAMIC=FALSE
    export MKL_NUM_THREADS=1
fi
module list

# RPATH integre au binaire depuis le CMakeLists.txt simplifie — pas besoin
# d'exporter LD_LIBRARY_PATH en temps normal. Filet de securite conserve
# pour le cas ou un binaire serait deplace hors de son build dir.
export LD_LIBRARY_PATH="$HOME/dependances_pour_serveur/install/suitesparse-install/lib:${LD_LIBRARY_PATH:-}"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-${OMP_NUM_THREADS:-1}}"
export OMP_PROC_BIND="${OMP_PROC_BIND:-close}"
export OMP_PLACES="${OMP_PLACES:-cores}"

MIN_BINARY="$MIN_BUILD_DIR/PhaseFieldFerro"

if [[ ! -x "$MIN_BINARY" ]]; then
    echo "[common_env] ATTENTION : binaire introuvable ou non exécutable : $MIN_BINARY" >&2
fi

export MIN_BUILD_DIR MIN_BINARY

echo "[common_env] variante=$MIN_VARIANT build_dir=$MIN_BUILD_DIR OMP_NUM_THREADS=$OMP_NUM_THREADS MKL=$MIN_USE_MKL"
