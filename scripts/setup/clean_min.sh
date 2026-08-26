#!/bin/bash
TARGET_DIR="$HOME/Min"
DRY_RUN=1
[[ "$1" == "--apply" ]] && DRY_RUN=0

function do_cmd() {
    if [ $DRY_RUN -eq 1 ]; then
        echo "[DRY-RUN] $@"
    else
        echo "[APPLY] $@"
        eval "$@"
    fi
}

echo "=== Réorganisation de $TARGET_DIR ==="

# 1. Création des dossiers
do_cmd "mkdir -p $TARGET_DIR/{scripts,logs,results,archives}"

# 2. (obsolète) Le renommage build->build_release / build-profiling->build_tracy
# n'a plus lieu d'être : les 4 dossiers de build (build, build_release,
# build_tracy, build_pardiso_tracy) sont maintenant nommés explicitement
# dès leur création par configure.sh — "build" seul désigne désormais le
# build PARDISO de prod, pas un ancien nom générique à renommer.

# 3. Rangement des scripts
do_cmd "mv $TARGET_DIR/*.sh $TARGET_DIR/*.py $TARGET_DIR/scripts/ 2>/dev/null"

# 4. Rangement des résultats
do_cmd "mv $TARGET_DIR/bench_* $TARGET_DIR/scaling.png $TARGET_DIR/results/ 2>/dev/null"
# (bench_scaling.sh a déjà été bougé par l'étape scripts, il ne reste que les dossiers bench_*)

# 5. Rangement des logs bruts SLURM
do_cmd "mv $TARGET_DIR/*.out $TARGET_DIR/run.sh.* $TARGET_DIR/logs/ 2>/dev/null"

# 6. Exfiltration des archives parasites
do_cmd "mv $TARGET_DIR/core/src/Physics.zip $TARGET_DIR/archives/ 2>/dev/null"

echo "Terminé. $([ $DRY_RUN -eq 1 ] && echo "Relancez avec --apply pour exécuter." || echo "Modifications appliquées.")"