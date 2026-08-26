#!/bin/bash
TARGET_DIR="$HOME/dependances_pour_serveur"
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

# 1. Création de la structure
do_cmd "mkdir -p $TARGET_DIR/{install,src,archive}"

# 2. Isolation des installations propres
do_cmd "mv $TARGET_DIR/eigen-install $TARGET_DIR/mimalloc-install $TARGET_DIR/suitesparse-install $TARGET_DIR/gmsh-sdk $TARGET_DIR/install/ 2>/dev/null"

# 3. Isolation de la source active de Tracy (version épinglée)
do_cmd "mv $TARGET_DIR/tracy-0.13.1-src $TARGET_DIR/src/ 2>/dev/null"

# 4. Mise aux archives de tout le reste (Zips, dossiers sources inutiles post-install, master incertain)
# On déplace d'abord les fichiers zips/tgz
do_cmd "mv $TARGET_DIR/*.zip $TARGET_DIR/*.tgz $TARGET_DIR/archive/ 2>/dev/null"
# Puis les dossiers sources qui ont déjà une version "install" ou qui sont obsolètes
do_cmd "mv $TARGET_DIR/eigen-3.4.0 $TARGET_DIR/mimalloc-master $TARGET_DIR/SuiteSparse-master $TARGET_DIR/tracy-master $TARGET_DIR/capstone-src $TARGET_DIR/archive/ 2>/dev/null"

echo "Terminé. $([ $DRY_RUN -eq 1 ] && echo "Relancez avec --apply pour exécuter." || echo "Modifications appliquées.")"