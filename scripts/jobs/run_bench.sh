#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=PhaseFieldFerro-bench
#SBATCH --output=%x-%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48

#SBATCH --time=12:00:00

# =========================================================================
# Usage :
#   sbatch run_bench.sh <config.toml> <out_dir>
# Exemples :
#   sbatch run_bench.sh ../input/config_bench_mechanics.toml bench_mechanics
#   sbatch run_bench.sh ../input/config_bench_full.toml      bench_full
# =========================================================================
CONFIG="${1:?Usage: sbatch run_bench.sh <config.toml> <out_dir>}"
OUT_DIR="${2:?Usage: sbatch run_bench.sh <config.toml> <out_dir>}"

# CORRECTION : $CONFIG et $OUT_DIR sont donnés relativement au dossier de
# soumission (~/Min), mais on fait ensuite `cd ~/Min/build_release` — sans résoudre
# les chemins en absolu avant, $CONFIG devient invalide (Datafile lève une
# exception dès le démarrage, échec quasi instantané et silencieux dans les
# logs "[TEMPS]"), et $OUT_DIR écrirait les résultats dans build/ au lieu
# du dossier attendu par l'utilisateur.
CONFIG="$(readlink -f "$CONFIG")"
if [[ ! -f "$CONFIG" ]]; then
    echo "[ERREUR] Fichier de config introuvable : $CONFIG" >&2
    exit 1
fi
mkdir -p "$OUT_DIR"
OUT_DIR="$(readlink -f "$OUT_DIR")"

# --- Environnement ---------------------------------------------------------
source /etc/profile.d/z99-modules.sh
module purge
module load mkl/latest
module list

# =========================================================================
# MKL : un seul runtime OpenMP dans le process (evite le conflit
# libgomp / libiomp5), threading géré par vos propres #pragma omp.
# =========================================================================
export MKL_THREADING_LAYER=GNU
export MKL_INTERFACE_LAYER=LP64
export MKL_DYNAMIC=FALSE
export MKL_NUM_THREADS=1

# CORRECTION (historique) : libcholmod.so vivait dans un prefix non-standard
# et n'etait pas trouve par le job SLURM sur w95/w96. Depuis le CMakeLists.txt
# simplifie, le RPATH est integre directement dans le binaire (BUILD_RPATH/
# INSTALL_RPATH) donc cette ligne n'est plus strictement necessaire — gardee
# en filet de securite si jamais le binaire est deplace hors de son build dir.
export LD_LIBRARY_PATH="$HOME/dependances_pour_serveur/install/suitesparse-install/lib:${LD_LIBRARY_PATH:-}"

# NOTE : OMP_NUM_THREADS / OMP_PROC_BIND / OMP_PLACES ne sont PAS fixés ici
# volontairement — bench_scaling.sh les définit lui-même à chaque run
# (variation contrôlée du nombre de threads, OMP_PROC_BIND=close), toute
# valeur fixée ici serait de toute façon écrasée. Les fixer ici porterait
# à confusion sur qui contrôle réellement le placement.

cd ~/Min/build_release

# CORRECTION : dos2unix défensif — bench_scaling.sh a déjà été retransféré/
# réédité avec des CRLF plusieurs fois (probablement lié à l'outil de
# transfert ou d'édition côté poste local). On le nettoie systématiquement
# ici plutôt que de compter sur le fait qu'il soit propre à chaque fois.
sed -i 's/\r$//' ~/Min/scripts/bench_scaling.sh
chmod +x ~/Min/scripts/bench_scaling.sh

echo "Job SLURM $SLURM_JOB_ID sur $(hostname), $SLURM_CPUS_PER_TASK coeurs alloués (exclusif)."
echo "Config: $CONFIG"
echo "Sortie: $OUT_DIR"

OUT_DIR="$OUT_DIR" ~/Min/scripts/bench_scaling.sh ./PhaseFieldFerro "$CONFIG"
bench_exit=$?

# CORRECTION : avant on affichait "Balayage terminé" même si bench_scaling.sh
# avait planté immédiatement (ex: CRLF dans son shebang) — le job SLURM
# finissait "COMPLETED" en 0 seconde sans que rien n'ait vraiment tourné.
# On vérifie maintenant explicitement le code de retour et le contenu du CSV.
if [[ $bench_exit -ne 0 ]]; then
    echo "[ERREUR] bench_scaling.sh a retourné le code $bench_exit — le balayage n'a probablement pas tourné." >&2
    exit 1
fi

n_data_lines=$(( $(wc -l < "$OUT_DIR/results.csv" 2>/dev/null || echo 0) - 1 ))
if [[ ! -f "$OUT_DIR/results.csv" || $n_data_lines -le 0 ]]; then
    echo "[ERREUR] $OUT_DIR/results.csv est vide ou absent — le balayage n'a produit aucune donnée." >&2
    exit 1
fi

echo "Balayage terminé avec succès : $n_data_lines lignes de données dans $OUT_DIR/results.csv"