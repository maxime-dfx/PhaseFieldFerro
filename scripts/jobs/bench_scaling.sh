#!/usr/bin/env bash
# =============================================================================
# bench_scaling.sh — Phase 0 : mesure de scalabilité avant refactorisation
#
# Fait varier OMP_NUM_THREADS et la politique NUMA, lance le binaire de
# simulation pour chaque combinaison, extrait le temps total ("[TEMPS]
# SimulationManager completed in : X ms") et écrit un CSV exploitable.
#
# Prérequis :
#   - Compiler le binaire normalement (le script ne compile rien).
#   - Dans le config.toml utilisé pour le benchmark : chrono_run = true
#     (sous la section [chrono]) — sinon la ligne [TEMPS] n'apparaît pas.
#   - Choisir un cas de taille "moyenne" : assez gros pour que le calcul
#     domine le bruit de démarrage, assez petit pour boucler sur ~10
#     configurations de threads en un temps raisonnable. Un run de
#     quelques dizaines de secondes à 1 thread est un bon ordre de grandeur.
#
# Usage :
#   ./bench_scaling.sh /chemin/vers/binaire /chemin/vers/config.toml
#
# Variables d'environnement optionnelles :
#   THREAD_LIST   liste des nombres de threads à tester (défaut ci-dessous)
#   REPEATS       nombre de répétitions par configuration (défaut 3)
#   OUT_DIR       dossier de sortie (défaut ./bench_results_<date>)
#   NUMA_MODES    liste de politiques numactl à tester (défaut ci-dessous)
# =============================================================================
set -euo pipefail

BIN="${1:?Usage: $0 <binaire> <config.toml>}"
CONFIG="${2:?Usage: $0 <binaire> <config.toml>}"

# --- Paramétrage ------------------------------------------------------------
# 24 = coeurs physiques (12/socket x 2), 48 = logiques (HT). On teste les deux
# côtés de cette frontière plus finement, c'est là que le comportement change.
THREAD_LIST="${THREAD_LIST:-1 2 4 8 12 16 20 24 28 32 40 48}"
REPEATS="${REPEATS:-3}"
OUT_DIR="${OUT_DIR:-./bench_results_$(date +%Y%m%d_%H%M%S)}"

# Politiques NUMA à comparer :
#   default    : aucune contrainte (comportement actuel du code)
#   local1sock : tout sur le socket 0 (cpu + mémoire) — référence "sans NUMA"
#   interleave : mémoire entrelacée sur les 2 noeuds — utile pour juger
#                l'impact du first-touch actuel avant de le corriger
NUMA_MODES="${NUMA_MODES:-default local1sock interleave}"

mkdir -p "$OUT_DIR"
CSV="$OUT_DIR/results.csv"
echo "numa_mode,threads,repeat,wall_ms,exit_code" > "$CSV"

if ! command -v numactl >/dev/null 2>&1; then
    echo "[ATTENTION] numactl introuvable — les modes 'local1sock' et 'interleave' seront sautés." >&2
    NUMA_MODES="default"
fi

echo "Binaire      : $BIN"
echo "Config       : $CONFIG"
echo "Threads      : $THREAD_LIST"
echo "NUMA modes   : $NUMA_MODES"
echo "Répétitions  : $REPEATS"
echo "Sortie       : $OUT_DIR"
echo

run_one() {
    local numa_mode="$1" threads="$2" rep="$3"
    local log_file="$OUT_DIR/log_${numa_mode}_t${threads}_r${rep}.txt"
    local prefix=()

    case "$numa_mode" in
        default)    prefix=() ;;
        local1sock) prefix=(numactl --cpunodebind=0 --membind=0) ;;
        interleave) prefix=(numactl --interleave=all) ;;
        *) echo "Mode NUMA inconnu: $numa_mode" >&2; return 1 ;;
    esac

    export OMP_NUM_THREADS="$threads"
    export OMP_PROC_BIND=close
    export OMP_PLACES=cores

    local exit_code=0
    "${prefix[@]}" "$BIN" "$CONFIG" > "$log_file" 2>&1 || exit_code=$?

    # Ligne attendue (avec codes couleur ANSI à retirer) :
    # [TEMPS] SimulationManager completed in  : 12345.6 ms
    local wall_ms
    wall_ms=$(sed -r 's/\x1B\[[0-9;]*[a-zA-Z]//g' "$log_file" \
        | grep -oP '(?<=SimulationManager completed in )[^:]*:\s*\K[0-9]+\.?[0-9]*' \
        | tail -1 || true)

    if [[ -z "$wall_ms" ]]; then
        echo "  [!] threads=$threads numa=$numa_mode rep=$rep : pas de temps trouvé (voir $log_file)" >&2
        wall_ms="NA"
    else
        echo "  threads=$threads numa=$numa_mode rep=$rep : ${wall_ms} ms"
    fi

    echo "${numa_mode},${threads},${rep},${wall_ms},${exit_code}" >> "$CSV"
}

for numa_mode in $NUMA_MODES; do
    echo "=== Mode NUMA: $numa_mode ==="
    for threads in $THREAD_LIST; do
        for rep in $(seq 1 "$REPEATS"); do
            run_one "$numa_mode" "$threads" "$rep"
        done
    done
    echo
done

echo "Terminé. Résultats bruts : $CSV"
echo "Logs individuels dans   : $OUT_DIR/"
echo
echo "Pour visualiser : python3 plot_scaling.py $CSV"