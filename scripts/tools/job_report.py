#!/usr/bin/env python3
"""
job_report.py — rapport d'efficacité CPU pour un ou plusieurs jobs SLURM.

Interroge `sacct` pour comparer le temps CPU réellement consommé au temps
CPU alloué (cpus-per-task x wall time), et signale les jobs qui gaspillent
des cores (utile pour repérer un mauvais OMP_NUM_THREADS, un solveur qui
reste mono-thread malgré l'allocation, une section I/O-bound, etc.)

Usage :
    python3 job_report.py <job_id> [<job_id> ...]
    python3 job_report.py --since 2026-08-01               # tous les jobs depuis une date
    python3 job_report.py --user $USER --since 2026-08-01
    python3 job_report.py 12345 --csv report.csv

Colonnes clés interprétées :
    Elapsed     temps mur du job
    NCPUS       cœurs alloués
    TotalCPU    temps CPU cumulé réellement utilisé (tous cœurs confondus)
    MaxRSS      pic mémoire résidente

Efficacité CPU = TotalCPU / (Elapsed * NCPUS)
    ~1.0  -> tous les cœurs alloués ont été utilisés en continu (bien)
    <0.5  -> plus de la moitié des cœurs alloués étaient inactifs en moyenne
             (sur-allocation, section séquentielle dominante, attente I/O...)
"""

from __future__ import annotations

import argparse
import csv as csv_module
import subprocess
import sys
from dataclasses import dataclass


def parse_slurm_time_to_seconds(t: str) -> float:
    """Parse un temps SLURM du type '1-02:03:04', '02:03:04' ou '03:04'."""
    if not t or t in ("", "Unknown", "INVALID"):
        return 0.0
    days = 0
    if "-" in t:
        days_str, t = t.split("-", 1)
        days = int(days_str)
    parts = [float(p) for p in t.split(":")]
    while len(parts) < 3:
        parts.insert(0, 0.0)
    h, m, s = parts[-3], parts[-2], parts[-1]
    return days * 86400 + h * 3600 + m * 60 + s


def parse_mem_to_mb(mem: str) -> float:
    """Parse une mémoire SLURM du type '1234K', '5.6M', '2G'."""
    if not mem:
        return 0.0
    mem = mem.strip()
    try:
        if mem.endswith("K"):
            return float(mem[:-1]) / 1024
        if mem.endswith("M"):
            return float(mem[:-1])
        if mem.endswith("G"):
            return float(mem[:-1]) * 1024
        return float(mem) / (1024 * 1024)  # bytes bruts, cas rare
    except ValueError:
        return 0.0


@dataclass
class JobRecord:
    job_id: str
    job_name: str
    state: str
    elapsed_s: float
    ncpus: int
    total_cpu_s: float
    max_rss_mb: float

    @property
    def cpu_efficiency(self) -> float:
        allocated = self.elapsed_s * self.ncpus
        if allocated <= 0:
            return 0.0
        return self.total_cpu_s / allocated


def fetch_jobs(job_ids: list[str] | None, user: str | None, since: str | None) -> list[JobRecord]:
    fields = "JobID,JobName,State,Elapsed,NCPUS,TotalCPU,MaxRSS"
    cmd = ["sacct", "--noheader", "--parsable2", f"--format={fields}"]

    if job_ids:
        cmd += ["-j", ",".join(job_ids)]
    if user:
        cmd += ["-u", user]
    if since:
        cmd += ["-S", since]
    if not job_ids and not user and not since:
        print("[job_report] Ni job_id, ni --user, ni --since fourni — utilisation "
              "de l'utilisateur courant sans filtre de date (peut être lent).",
              file=sys.stderr)

    try:
        proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
    except FileNotFoundError:
        print("[ERREUR] sacct introuvable — ce script doit tourner sur le cluster SLURM.",
              file=sys.stderr)
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"[ERREUR] sacct a échoué : {e.stderr}", file=sys.stderr)
        sys.exit(1)

    records: list[JobRecord] = []
    for line in proc.stdout.splitlines():
        parts = line.split("|")
        if len(parts) != 7:
            continue
        job_id, job_name, state, elapsed, ncpus, total_cpu, max_rss = parts

        # sacct renvoie une ligne par "step" (JobID, JobID.batch, JobID.extern...).
        # On ne garde que la ligne principale (JobID pur, sans suffixe .xxx),
        # sauf pour TotalCPU/MaxRSS où le step .batch est souvent plus fiable
        # -- on prend simplement le max observé par job_id racine.
        root_id = job_id.split(".")[0]
        is_root = (job_id == root_id)

        elapsed_s = parse_slurm_time_to_seconds(elapsed)
        total_cpu_s = parse_slurm_time_to_seconds(total_cpu)
        max_rss_mb = parse_mem_to_mb(max_rss)

        existing = next((r for r in records if r.job_id == root_id), None)
        if existing is None:
            try:
                ncpus_int = int(ncpus)
            except ValueError:
                ncpus_int = 0
            records.append(JobRecord(
                job_id=root_id, job_name=job_name, state=state,
                elapsed_s=elapsed_s, ncpus=ncpus_int,
                total_cpu_s=total_cpu_s, max_rss_mb=max_rss_mb,
            ))
        else:
            # Complète avec les infos des sous-steps si la ligne racine
            # ne les avait pas (cas fréquent : TotalCPU/MaxRSS vides sur
            # la ligne racine, présents sur .batch)
            existing.total_cpu_s = max(existing.total_cpu_s, total_cpu_s)
            existing.max_rss_mb = max(existing.max_rss_mb, max_rss_mb)
            if is_root:
                existing.state = state
                existing.job_name = job_name
                existing.elapsed_s = elapsed_s or existing.elapsed_s

    return records


def print_report(records: list[JobRecord]) -> None:
    if not records:
        print("Aucun job trouvé.")
        return

    header = f"{'JobID':>10} {'Name':<28} {'State':<12} {'Elapsed':>10} {'NCPUS':>6} {'TotalCPU':>10} {'Eff.':>6} {'MaxRSS(MB)':>11}"
    print(header)
    print("-" * len(header))
    for r in sorted(records, key=lambda r: r.job_id):
        eff = r.cpu_efficiency
        flag = ""
        if r.state == "COMPLETED" and r.ncpus > 1:
            if eff < 0.3:
                flag = "  <-- très faible efficacité, vérifier le parallélisme"
            elif eff < 0.6:
                flag = "  <-- efficacité moyenne"
        print(f"{r.job_id:>10} {r.job_name[:28]:<28} {r.state:<12} "
              f"{r.elapsed_s/60:>9.1f}m {r.ncpus:>6} {r.total_cpu_s/60:>9.1f}m "
              f"{eff:>6.1%} {r.max_rss_mb:>11.0f}{flag}")


def write_csv(records: list[JobRecord], path: str) -> None:
    with open(path, "w", newline="") as f:
        writer = csv_module.writer(f)
        writer.writerow(["job_id", "job_name", "state", "elapsed_s", "ncpus",
                          "total_cpu_s", "cpu_efficiency", "max_rss_mb"])
        for r in sorted(records, key=lambda r: r.job_id):
            writer.writerow([r.job_id, r.job_name, r.state, r.elapsed_s, r.ncpus,
                              r.total_cpu_s, f"{r.cpu_efficiency:.4f}", r.max_rss_mb])
    print(f"\nCSV écrit : {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("job_ids", nargs="*", help="ID(s) de job SLURM à inspecter")
    parser.add_argument("--user", type=str, default=None,
                         help="Filtrer par utilisateur (sacct -u)")
    parser.add_argument("--since", type=str, default=None,
                         help="Filtrer depuis une date (sacct -S), ex: 2026-08-01")
    parser.add_argument("--csv", type=str, default=None,
                         help="Écrire le rapport dans un CSV en plus de l'affichage")
    args = parser.parse_args()

    records = fetch_jobs(args.job_ids or None, args.user, args.since)
    print_report(records)
    if args.csv:
        write_csv(records, args.csv)


if __name__ == "__main__":
    main()
