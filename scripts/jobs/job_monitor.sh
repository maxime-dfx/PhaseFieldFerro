#!/usr/bin/env bash
# =============================================================================
# job_monitor.sh — tableau de bord live pour un job SLURM en cours.
#
# Affiche, rafraîchi toutes les N secondes :
#   - l'état du job (squeue : état, temps écoulé, nœud)
#   - l'usage CPU/mémoire réel sur le nœud (via ssh + ps, si le job tourne)
#   - les dernières lignes du log SLURM (%x-%j.out)
#
# Usage :
#   ./job_monitor.sh <job_id> [--log <chemin>] [--interval <secondes>]
#
# Exemple :
#   ./job_monitor.sh 12345
#   ./job_monitor.sh 12345 --log ~/Min/logs/PhaseFieldFerro-12345.out --interval 10
# =============================================================================
set -euo pipefail

JOB_ID="${1:?Usage: $0 <job_id> [--log <chemin>] [--interval <secondes>]}"
shift || true

LOG_PATH=""
INTERVAL=5

while [[ $# -gt 0 ]]; do
    case "$1" in
        --log) LOG_PATH="$2"; shift 2 ;;
        --interval) INTERVAL="$2"; shift 2 ;;
        *) echo "Argument inconnu: $1" >&2; exit 1 ;;
    esac
done

if ! command -v squeue >/dev/null 2>&1; then
    echo "[ERREUR] squeue introuvable — ce script doit tourner sur le cluster SLURM." >&2
    exit 1
fi

# Si le chemin du log n'est pas donné, on tente de le deviner depuis le
# format par défaut des scripts run_*.sh (%x-%j.out dans le dossier courant
# ou ~/Min).
guess_log_path() {
    local job_name
    job_name=$(squeue -h -j "$JOB_ID" -o "%j" 2>/dev/null || true)
    [[ -z "$job_name" ]] && return 1
    for dir in . "$HOME/Min" "$HOME/Min/logs"; do
        local candidate="$dir/${job_name}-${JOB_ID}.out"
        [[ -f "$candidate" ]] && { echo "$candidate"; return 0; }
    done
    return 1
}

if [[ -z "$LOG_PATH" ]]; then
    LOG_PATH="$(guess_log_path || true)"
fi

clear_screen() { printf '\033[2J\033[H'; }

while true; do
    clear_screen
    echo "=== job_monitor.sh — job $JOB_ID — $(date '+%Y-%m-%d %H:%M:%S') ==="
    echo

    # --- État SLURM -----------------------------------------------------
    STATE_LINE=$(squeue -h -j "$JOB_ID" -o "%T|%M|%L|%N|%C" 2>/dev/null || true)
    if [[ -z "$STATE_LINE" ]]; then
        echo "[squeue] Job $JOB_ID introuvable dans la file — probablement terminé."
        echo "         Vérifiez avec : sacct -j $JOB_ID --format=JobID,State,Elapsed,MaxRSS,TotalCPU"
        break
    fi
    IFS='|' read -r STATE ELAPSED REMAIN NODE NCPUS <<< "$STATE_LINE"
    echo "État        : $STATE"
    echo "Temps écoulé: $ELAPSED   (restant estimé: $REMAIN)"
    echo "Nœud        : $NODE"
    echo "CPUs alloués: $NCPUS"
    echo

    # --- Usage réel sur le nœud (si RUNNING) -----------------------------
    if [[ "$STATE" == "RUNNING" && -n "$NODE" && "$NODE" != "(null)" ]]; then
        echo "--- Usage CPU/mémoire sur $NODE (top 5 process du job) ---"
        # sstat donne l'usage agrégé du job directement, sans ssh.
        sstat -j "${JOB_ID}.batch" \
            --format=JobID,AveCPU,MaxRSS,MaxVMSize,NTasks 2>/dev/null \
            || echo "  (sstat indisponible ou job pas encore mesurable)"
        echo
    fi

    # --- Log applicatif ----------------------------------------------------
    if [[ -n "$LOG_PATH" && -f "$LOG_PATH" ]]; then
        echo "--- Dernières lignes de $LOG_PATH ---"
        tail -n 15 "$LOG_PATH"
    else
        echo "(log non trouvé — précisez avec --log <chemin>)"
    fi

    echo
    echo "[rafraîchi toutes les ${INTERVAL}s — Ctrl+C pour quitter]"
    sleep "$INTERVAL"
done
