#!/usr/bin/env python3
"""
plot_scaling.py — visualise les résultats produits par bench_scaling.sh

Usage:
    python3 plot_scaling.py bench_results_XXXXXXXX/results.csv [--out scaling.png]

Produit un PNG avec :
  - temps mur médian vs nombre de threads (par mode NUMA)
  - speedup relatif au run 1-thread du même mode, vs speedup idéal (linéaire)
  - efficacité parallèle (speedup / threads), pour repérer où ça décroche

Lignes verticales à 24 (coeurs physiques) et 48 (logiques, HT) pour situer
visuellement les seuils "coeurs physiques épuisés" / "HT épuisé".
"""
import sys
import argparse
import csv
from collections import defaultdict
import statistics as stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_results(path):
    # data[numa_mode][threads] = [wall_ms, ...]
    data = defaultdict(lambda: defaultdict(list))
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["wall_ms"] in ("NA", "", None):
                continue
            if row["exit_code"] != "0":
                # on garde quand même mais on prévient
                print(f"[warn] exit_code={row['exit_code']} pour "
                      f"numa={row['numa_mode']} threads={row['threads']} "
                      f"repeat={row['repeat']} — vérifier le log", file=sys.stderr)
            data[row["numa_mode"]][int(row["threads"])].append(float(row["wall_ms"]))
    return data


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv_path")
    ap.add_argument("--out", default="scaling.png")
    ap.add_argument("--physical-cores", type=int, default=24)
    ap.add_argument("--logical-cores", type=int, default=48)
    args = ap.parse_args()

    data = load_results(args.csv_path)
    if not data:
        print("Aucune donnée exploitable dans le CSV.", file=sys.stderr)
        sys.exit(1)

    fig, (ax_time, ax_speedup, ax_eff) = plt.subplots(1, 3, figsize=(18, 5.5))

    for numa_mode, per_threads in sorted(data.items()):
        threads_sorted = sorted(per_threads.keys())
        medians = [stats.median(per_threads[t]) for t in threads_sorted]

        # Référence pour le speedup : le plus petit nombre de threads mesuré
        # pour ce mode (normalement 1, mais on s'adapte si absent).
        t0 = threads_sorted[0]
        base = stats.median(per_threads[t0])

        speedups = [base / m for m in medians]
        efficiencies = [s / (t / t0) for s, t in zip(speedups, threads_sorted)]

        ax_time.plot(threads_sorted, medians, marker="o", label=numa_mode)
        ax_speedup.plot(threads_sorted, speedups, marker="o", label=numa_mode)
        ax_eff.plot(threads_sorted, efficiencies, marker="o", label=numa_mode)

    # Speedup idéal (linéaire) pris sur l'axe des threads du 1er mode tracé
    any_mode = next(iter(data.values()))
    all_threads = sorted(any_mode.keys())
    ideal = [t / all_threads[0] for t in all_threads]
    ax_speedup.plot(all_threads, ideal, "k--", alpha=0.5, label="idéal (linéaire)")
    ax_eff.axhline(1.0, color="k", linestyle="--", alpha=0.5, label="idéal")

    for ax in (ax_time, ax_speedup, ax_eff):
        ax.axvline(args.physical_cores, color="gray", linestyle=":", alpha=0.7)
        ax.axvline(args.logical_cores, color="gray", linestyle=":", alpha=0.7)
        ax.set_xlabel("threads (OMP_NUM_THREADS)")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)

    ax_time.set_ylabel("temps mur médian (ms)")
    ax_time.set_title("Temps d'exécution")
    ax_time.set_yscale("log")

    ax_speedup.set_ylabel("speedup")
    ax_speedup.set_title("Speedup vs 1 thread")

    ax_eff.set_ylabel("efficacité (speedup / threads)")
    ax_eff.set_title("Efficacité parallèle")
    ax_eff.set_ylim(0, 1.15)

    fig.suptitle(
        f"Scalabilité — pointillés gris = {args.physical_cores} coeurs physiques "
        f"/ {args.logical_cores} logiques (HT)"
    )
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"Graphique écrit : {args.out}")


if __name__ == "__main__":
    main()