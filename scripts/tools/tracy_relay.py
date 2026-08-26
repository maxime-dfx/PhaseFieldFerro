#!/usr/bin/env python3
"""
tracy_relay.py — relaie un port local vers le port Tracy d'un job SLURM
en cours, via `srun --overlap`, pour profiler en direct un job qui tourne
sur un nœud de calcul sans accès réseau direct.

Corrections par rapport à la version précédente (test.py) :
  - le job_id n'est plus codé en dur : passé en argument (--job-id),
    obligatoire, ce qui évitait un relais silencieux vers le mauvais job
    si le script était réutilisé sans y repenser.
  - ports configurables en ligne de commande.
  - le process srun est vérifié au démarrage (erreur claire si le job_id
    n'existe pas / n'est pas RUNNING, plutôt qu'un relais qui accepte des
    connexions mais ne transmet jamais rien).

Usage :
    python3 tracy_relay.py --job-id 781
    python3 tracy_relay.py --job-id 781 --local-port 8086 --remote-port 8086
"""
import argparse
import socket
import subprocess
import sys
import threading


def check_job_running(job_id: str) -> None:
    try:
        proc = subprocess.run(
            ["squeue", "-h", "-j", job_id, "-o", "%T"],
            capture_output=True, text=True, timeout=10,
        )
    except FileNotFoundError:
        print("[tracy_relay] squeue introuvable — ce script doit tourner "
              "sur le cluster SLURM.", file=sys.stderr)
        sys.exit(1)

    state = proc.stdout.strip()
    if not state:
        print(f"[tracy_relay] job {job_id} introuvable dans la file "
              f"(terminé ou id invalide).", file=sys.stderr)
        sys.exit(1)
    if state != "RUNNING":
        print(f"[tracy_relay] ATTENTION : job {job_id} est dans l'état "
              f"'{state}' (pas RUNNING) — le relais risque de ne rien "
              f"transmettre tant que le job n'est pas actif.", file=sys.stderr)


def handle(conn: socket.socket, job_id: str, remote_port: int) -> None:
    remote_snippet = (
        f"import socket,sys,shutil,threading\n"
        f"s=socket.create_connection(('localhost',{remote_port}))\n"
        f"def pipe(a,b):\n"
        f"    while True:\n"
        f"        d=a.read(1)\n"
        f"        if not d: break\n"
        f"        b.write(d); b.flush()\n"
        f"t=threading.Thread(target=lambda: shutil.copyfileobj(s.makefile('rb'), sys.stdout.buffer))\n"
        f"t.start()\n"
        f"shutil.copyfileobj(sys.stdin.buffer, s.makefile('wb'))\n"
    )
    proc = subprocess.Popen(
        ["srun", f"--jobid={job_id}", "--overlap", "python3", "-c", remote_snippet],
        stdin=subprocess.PIPE, stdout=subprocess.PIPE,
    )

    def forward(src: socket.socket, dst) -> None:
        try:
            while True:
                d = src.recv(4096)
                if not d:
                    break
                dst.write(d)
                dst.flush()
        finally:
            dst.close()

    threading.Thread(target=forward, args=(conn, proc.stdin), daemon=True).start()
    while True:
        d = proc.stdout.read(4096)
        if not d:
            break
        conn.sendall(d)
    conn.close()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--job-id", required=True, help="ID du job SLURM ciblé (obligatoire)")
    parser.add_argument("--local-port", type=int, default=8086)
    parser.add_argument("--remote-port", type=int, default=8086)
    args = parser.parse_args()

    check_job_running(args.job_id)

    srv = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    srv.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    srv.bind(("0.0.0.0", args.local_port))
    srv.listen(5)
    print(f"[tracy_relay] relais actif sur :{args.local_port} -> "
          f"job {args.job_id} port {args.remote_port} (Ctrl+C pour arrêter)")
    try:
        while True:
            conn, _ = srv.accept()
            threading.Thread(target=handle, args=(conn, args.job_id, args.remote_port),
                              daemon=True).start()
    except KeyboardInterrupt:
        print("\n[tracy_relay] arrêt.")


if __name__ == "__main__":
    main()
