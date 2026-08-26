#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=PhaseFieldFerro_scaling
#SBATCH --output=%x-%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=24:00:00

# ---------------------------------------------------------------------
# Benchmark de scalabilité : on alloue 64 cores UNE SEULE FOIS, puis on
# lance le binaire plusieurs fois en variant OMP_NUM_THREADS.
# Avantage : pas de temps d'attente en queue entre chaque run.
# Attention : le run à 64 threads peut être influencé par le fait qu'on
# a déjà "chauffé" le nœud (cache, fréquence CPU) — répète si besoin.
# ---------------------------------------------------------------------

source /etc/profile.d/z99-modules.sh
module purge
module load mkl/latest
module list

export OMP_PROC_BIND=close
export OMP_PLACES=cores
export MKL_THREADING_LAYER=GNU
export MKL_INTERFACE_LAYER=LP64
export MKL_DYNAMIC=FALSE
export MKL_NUM_THREADS=1

cd ~/Min/build_pardiso_tracy

BINARY=./PhaseFieldFerro
CONFIG=../input/config.toml
REPEATS=3
CORE_LIST="1 2 4 8 16 24 32 48 64"

RESULTS_CSV="scaling_results_${SLURM_JOB_ID}.csv"
echo "cores,repeat,time_seconds" > "$RESULTS_CSV"

echo "Job SLURM $SLURM_JOB_ID sur $(hostname), $SLURM_CPUS_PER_TASK coeurs alloués."
echo "Liste des cores testés : $CORE_LIST"
echo ""

for n in $CORE_LIST; do
    if [ "$n" -gt "$SLURM_CPUS_PER_TASK" ]; then
        echo "⚠️  Skip cores=$n (dépasse les $SLURM_CPUS_PER_TASK cores alloués)"
        continue
    fi

    export OMP_NUM_THREADS=$n

    for r in $(seq 1 $REPEATS); do
        echo ">>> cores=$n  run $r/$REPEATS"
        start=$(date +%s.%N)

        "$BINARY" "$CONFIG" > "run_cores${n}_rep${r}.log" 2>&1

        end=$(date +%s.%N)
        elapsed=$(echo "$end - $start" | bc)
        echo "    temps = ${elapsed}s"
        echo "$n,$r,$elapsed" >> "$RESULTS_CSV"
    done
done

echo ""
echo "Terminé. Résultats dans $RESULTS_CSV"