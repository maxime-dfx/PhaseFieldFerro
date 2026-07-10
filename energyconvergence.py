#!/usr/bin/env python3
"""
check_energy_convergence.py
============================
Lit le fichier energies.csv produit par Diagnostics::append_csv (voir
Simulation::extract_and_save_results) et verifie deux proprietes de base
attendues d'un solveur de gradient-flow correct :

  1. L'energie totale doit decroitre de facon MONOTONE au cours du temps
     (propriete fondamentale d'un flot de gradient dissipatif).
  2. Aucune valeur NaN/Inf, aucun saut brutal (signe d'une instabilite
     numerique ou d'un bug de signe comme celui trouve dans
     PolarizationAssembler).

Usage
-----
    python3 check_energy_convergence.py /chemin/vers/energies.csv

Si aucune colonne "total" n'existe dans le CSV, le script additionne
automatiquement toutes les colonnes numeriques dont le nom contient
"energy"/"energie" (insensible a la casse) pour reconstruire un total
approximatif -- adapte les noms de colonnes ci-dessous si besoin.
"""

import sys
import csv
import numpy as np

ENERGY_KEYWORDS = ["energy", "energie", "u_", "w_", "chi", "elec", "surf", "gradient"]

def load_csv(path):
    with open(path, "r") as f:
        reader = csv.reader(f)
        header = next(reader)
        rows = [row for row in reader if row]
    header = [h.strip() for h in header]
    data = np.array(rows, dtype=float)
    return header, data

def find_time_column(header):
    for i, h in enumerate(header):
        if h.lower() in ("time", "t", "temps"):
            return i
    return 0  # fallback : premiere colonne

def find_total_energy_column(header):
    for i, h in enumerate(header):
        if h.lower() in ("total", "total_energy", "energie_totale", "h_total"):
            return i
    return None

def reconstruct_total(header, data):
    cols = [i for i, h in enumerate(header)
            if any(k in h.lower() for k in ENERGY_KEYWORDS)]
    if not cols:
        return None, []
    total = data[:, cols].sum(axis=1)
    names = [header[i] for i in cols]
    return total, names

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 check_energy_convergence.py energies.csv")
        sys.exit(1)

    path = sys.argv[1]
    header, data = load_csv(path)
    print(f"Colonnes detectees : {header}")
    print(f"Nombre de lignes (pas de temps enregistres) : {data.shape[0]}")
    print("-" * 70)

    t_idx = find_time_column(header)
    time = data[:, t_idx]

    total_idx = find_total_energy_column(header)
    if total_idx is not None:
        total_energy = data[:, total_idx]
        source = f"colonne '{header[total_idx]}'"
    else:
        total_energy, used_cols = reconstruct_total(header, data)
        if total_energy is None:
            print("ERREUR : impossible de trouver ou reconstruire une energie totale.")
            print("Adapte ENERGY_KEYWORDS ou fournis directement le bon indice de colonne.")
            sys.exit(1)
        source = f"somme des colonnes {used_cols}"

    print(f"Energie totale prise depuis : {source}")

    # ---- Test 1 : NaN / Inf ----
    n_bad = np.sum(~np.isfinite(total_energy))
    if n_bad > 0:
        print(f"[ECHEC] {n_bad} valeurs non finies (NaN/Inf) detectees dans l'energie totale !")
    else:
        print("[OK] Aucune valeur NaN/Inf detectee.")

    # ---- Test 2 : monotonie ----
    diffs = np.diff(total_energy)
    n_increases = np.sum(diffs > 1e-10)  # tolerance numerique
    if n_increases == 0:
        print("[OK] Energie totale strictement decroissante (ou constante) -- comportement attendu.")
    else:
        worst_idx = np.argmax(diffs)
        print(f"[ATTENTION] {n_increases} augmentations d'energie detectees sur {len(diffs)} pas.")
        print(f"            Plus grande augmentation : +{diffs[worst_idx]:.6e}")
        print(f"            entre t={time[worst_idx]:.6g} et t={time[worst_idx+1]:.6g}")
        print("            -> signe possible d'un bug de signe dans l'assemblage")
        print("               (voir le cas PolarizationAssembler corrige) ou d'un dt trop grand.")

    # ---- Test 3 : saut brutal relatif ----
    rel_jumps = np.abs(diffs) / (np.abs(total_energy[:-1]) + 1e-12)
    n_jumps = np.sum(rel_jumps > 0.5)  # >50% de variation relative en un seul pas
    if n_jumps > 0:
        print(f"[ATTENTION] {n_jumps} sauts relatifs > 50% detectes -- verifie dt/adaptativite.")
    else:
        print("[OK] Pas de saut brutal (> 50% relatif) entre deux pas de temps consecutifs.")

    print("-" * 70)
    print(f"Energie initiale : {total_energy[0]:.6e}")
    print(f"Energie finale   : {total_energy[-1]:.6e}")
    print(f"Variation totale : {total_energy[-1] - total_energy[0]:.6e}")

if __name__ == "__main__":
    main()