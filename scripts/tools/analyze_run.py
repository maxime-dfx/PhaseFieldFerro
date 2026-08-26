#!/usr/bin/env python3
"""
analyze_run.py — analyse rapide d'un (ou plusieurs) energies.csv produit
par PhaseFieldFerro, pour vérifier la santé d'un run sans avoir à ouvrir
le CSV à la main.

Détecte :
  - valeurs NaN / infinies
  - sauts brusques d'un LoadStep à l'autre (> seuil, potentielle
    instabilité numérique ou divergence du solveur)
  - un run qui s'arrête plus tôt que prévu (moins de load steps que la
    référence, si --expect-steps est donné)
  - colonnes constantes suspectes (ex: Surface toujours à 0, signe que la
    fissure n'a jamais progressé)

Usage :
    python3 analyze_run.py path/to/energies.csv
    python3 analyze_run.py paper_experiments_output/*/run_*/energies.csv
    python3 analyze_run.py energies.csv --column Surface --jump-threshold 0.2
    python3 analyze_run.py energies.csv --expect-steps 200
"""

from __future__ import annotations

import argparse
import csv
import glob
import math
import sys
from pathlib import Path


def load_csv(path: Path) -> tuple[list[str], dict[str, list[float]]]:
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames or []
        columns: dict[str, list[float]] = {name: [] for name in fieldnames}
        for row in reader:
            for name in fieldnames:
                raw = row.get(name, "")
                try:
                    columns[name].append(float(raw))
                except (ValueError, TypeError):
                    columns[name].append(math.nan)
    return fieldnames, columns


def analyze_one(path: Path, jump_threshold: float, expect_steps: int | None,
                 watch_column: str | None) -> list[str]:
    """Retourne la liste des problèmes détectés (vide = tout va bien)."""
    issues: list[str] = []

    if not path.exists():
        return [f"fichier introuvable : {path}"]

    try:
        fieldnames, columns = load_csv(path)
    except Exception as e:
        return [f"impossible de lire le CSV ({e})"]

    if not fieldnames:
        return ["CSV vide ou sans en-tête"]

    n_rows = len(next(iter(columns.values()), []))
    if n_rows == 0:
        return ["aucune ligne de données"]

    if expect_steps is not None and n_rows < expect_steps:
        issues.append(f"seulement {n_rows} lignes, attendu >= {expect_steps} "
                       f"(run probablement interrompu prématurément)")

    for col_name, values in columns.items():
        nan_idx = next((i for i, v in enumerate(values) if math.isnan(v)), None)
        inf_idx = next((i for i, v in enumerate(values) if math.isinf(v)), None)
        nan_count = sum(1 for v in values if math.isnan(v))
        inf_count = sum(1 for v in values if math.isinf(v))
        if nan_count:
            issues.append(f"colonne '{col_name}' : {nan_count} valeur(s) NaN "
                           f"(1ère occurrence à la ligne {nan_idx + 1})")
        if inf_count:
            issues.append(f"colonne '{col_name}' : {inf_count} valeur(s) infinie(s) "
                           f"(1ère occurrence à la ligne {inf_idx + 1})")

    # Colonnes numériques "intéressantes" pour la détection de sauts : tout
    # sauf LoadStep lui-même (qui est monotone par construction).
    columns_to_watch = [watch_column] if watch_column else [
        c for c in fieldnames if c.lower() not in ("loadstep", "step", "iteration")
    ]

    for col_name in columns_to_watch:
        if col_name not in columns:
            issues.append(f"colonne demandée '{col_name}' absente du CSV "
                           f"(colonnes disponibles : {', '.join(fieldnames)})")
            continue
        values = columns[col_name]
        clean = [v for v in values if not (math.isnan(v) or math.isinf(v))]
        if len(clean) < 3:
            continue
        value_range = max(clean) - min(clean)
        if value_range == 0:
            pass  # géré plus bas par la détection de colonne constante
        else:
            # Un "saut" est jugé suspect s'il est à la fois :
            #  - une fraction significative de l'amplitude totale du run
            #    (jump_threshold), et
            #  - un net outlier par rapport aux autres variations pas-à-pas
            #    (pour ne pas signaler une progression régulière et attendue,
            #    où chaque pas est "grand" simplement parce que le run a peu
            #    de points ou démarre près de zéro).
            diffs = [abs(values[i] - values[i - 1]) for i in range(1, len(values))
                      if not (math.isnan(values[i]) or math.isnan(values[i - 1]))]
            diffs_sorted = sorted(diffs)
            median_diff = diffs_sorted[len(diffs_sorted) // 2] if diffs_sorted else 0.0
            outlier_floor = max(3 * median_diff, 1e-12)

            for i in range(1, len(values)):
                a, b = values[i - 1], values[i]
                if math.isnan(a) or math.isnan(b):
                    continue
                d = abs(b - a)
                rel = d / value_range
                if rel > jump_threshold and d > outlier_floor:
                    issues.append(f"colonne '{col_name}' : saut brusque de {a:.4g} à "
                                   f"{b:.4g} entre les lignes {i} et {i+1} "
                                   f"({rel:.1%} de l'amplitude du run, et "
                                   f"{d/median_diff:.1f}x le pas médian)")

        if all(v == clean[0] for v in clean):
            issues.append(f"colonne '{col_name}' : constante à {clean[0]:.4g} "
                           f"sur tout le run (rien ne semble évoluer)")

    return issues


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("paths", nargs="+",
                         help="Chemin(s) vers energies.csv (glob supporté)")
    parser.add_argument("--column", type=str, default=None,
                         help="N'analyser qu'une seule colonne (défaut : toutes sauf LoadStep)")
    parser.add_argument("--jump-threshold", type=float, default=0.3,
                         help="Seuil relatif de saut entre deux lignes considéré "
                              "comme suspect (défaut : 0.3 = 30%% de l'amplitude max)")
    parser.add_argument("--expect-steps", type=int, default=None,
                         help="Nombre minimum de lignes attendu (détecte un arrêt prématuré)")
    args = parser.parse_args()

    all_paths: list[Path] = []
    for p in args.paths:
        expanded = glob.glob(p, recursive=True)
        all_paths.extend(Path(e) for e in (expanded or [p]))

    if not all_paths:
        print("Aucun fichier trouvé.", file=sys.stderr)
        sys.exit(1)

    n_ok, n_issues = 0, 0
    for path in sorted(all_paths):
        issues = analyze_one(path, args.jump_threshold, args.expect_steps, args.column)
        if issues:
            n_issues += 1
            print(f"\n[!] {path}")
            for issue in issues:
                print(f"    - {issue}")
        else:
            n_ok += 1
            print(f"[OK] {path}")

    print(f"\n{n_ok} run(s) sans anomalie détectée, {n_issues} run(s) avec au moins un problème.")
    sys.exit(1 if n_issues else 0)


if __name__ == "__main__":
    main()
