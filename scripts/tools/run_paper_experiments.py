#!/usr/bin/env python3
import os
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

# --- Configuration des chemins ---
PROJECT_ROOT = Path(__file__).resolve().parent
BASE_CONFIG = PROJECT_ROOT / "input" / "config.toml"
OUTPUT_ROOT = PROJECT_ROOT / "paper_experiments_output"
EXECUTABLE = PROJECT_ROOT / "build_pardiso_petsc" / "PhaseFieldFerro"
SLURM_DIR = OUTPUT_ROOT / "slurm"

@dataclass
class Experiment:
    name: str
    fracture_mode: str
    E: float
    freeze_polarization: bool = False

    def get_voltage(self, Lx: float) -> float:
        # V = -E * L (Défini dans la Section 3.3 du papier)
        return -self.E * Lx

def build_experiment_list() -> list[Experiment]:
    """Définit toutes les expériences du papier d'Abdollahi & Arias (2011)."""
    exps = []
    
    # 1. Chargement mécanique pur (Section 3.2)
    exps.append(Experiment("sec3_2_1_mech_permeable", "PERMEABLE", 0.0))
    exps.append(Experiment("sec3_2_2_mech_impermeable", "IMPERMEABLE", 0.0))
    
    # 2. Modèles de référence (Single-phase, polarisation désactivée)
    exps.append(Experiment("single_phase_perm", "PERMEABLE", 0.0, True))
    exps.append(Experiment("single_phase_imp", "IMPERMEABLE", 0.0, True))
    
    # 3. Chargement électromécanique (Section 3.3)
    # Fissures perméables
    for e_val in [1e-3, 3e-3, 4e-3, 5e-3, -1e-3, -3e-3, -4e-3, -5e-3]:
        exps.append(Experiment(f"sec3_3_perm_{e_val}", "PERMEABLE", e_val))
        
    # Fissures imperméables
    for e_val in [1e-3, 2e-3, 4e-3, 5e-3, -1e-3, -2e-3, -4e-3, -5e-3]:
        exps.append(Experiment(f"sec3_3_imp_{e_val}", "IMPERMEABLE", e_val))
        
    return exps

def extract_lx(config_text: str) -> float:
    """Extrait la longueur du domaine L_x pour calculer le potentiel."""
    match = re.search(r'(?m)^L_x\s*=\s*([\-0-9.eE]+)', config_text)
    return float(match.group(1)) if match else 200.0

def generate_config(base_text: str, exp: Experiment, output_dir: Path) -> str:
    """Modifie les paramètres du config.toml de façon robuste."""
    text = base_text
    lx = extract_lx(text)
    
    # 1. Remplacement du dossier de sortie
    text = re.sub(r'(?m)^(OutputDir\s*=\s*")[^"]*(")', 
                  repl=lambda m: f'{m.group(1)}{output_dir.as_posix()}{m.group(2)}', string=text)
                  
    # 2. Remplacement du mode de fracture
    text = re.sub(r'(?m)^(mode\s*=\s*")[^"]*(")', 
                  repl=lambda m: f'{m.group(1)}{exp.fracture_mode}{m.group(2)}', string=text)
                  
    # 3. Activation/Désactivation de la polarisation
    pol_state = "false" if exp.freeze_polarization else "true"
    text = re.sub(r'(?m)^(enable_polarization\s*=\s*)\w+', 
                  repl=lambda m: f'{m.group(1)}{pol_state}', string=text)
                  
    # 4. Modification ROBUSTE de la tension V sur le bord droit
    voltage_str = f"{exp.get_voltage(lx):.10g}"
    
    # On découpe le fichier TOML en blocs à chaque [[boundary_rule]]
    blocks = text.split("[[boundary_rule]]")
    
    # On cherche le bloc qui gère "phi" sur le bord "right"
    for i in range(1, len(blocks)):
        if 'field = "phi"' in blocks[i] and 'edge_name = "right"' in blocks[i]:
            # On remplace uniquement la ligne "val = ..." dans CE bloc spécifique
            blocks[i] = re.sub(r'(?m)^(val\s*=\s*)([\-0-9.eE]+)', 
                               repl=lambda m: f"{m.group(1)}{voltage_str}", string=blocks[i])
            break # On a trouvé et modifié le bon bloc, on s'arrête
            
    # On recompose le fichier complet
    text = "[[boundary_rule]]".join(blocks)
    
    return text

def setup_and_submit_slurm():
    base_text = BASE_CONFIG.read_text(encoding="utf-8")
    experiments = build_experiment_list()
    
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    SLURM_DIR.mkdir(parents=True, exist_ok=True)
    
    exp_dirs = []
    
    # 1. Génération de tous les fichiers d'entrée
    for exp in experiments:
        exp_dir = OUTPUT_ROOT / exp.name
        exp_dir.mkdir(parents=True, exist_ok=True)
        
        new_config = generate_config(base_text, exp, exp_dir)
        (exp_dir / "config.toml").write_text(new_config, encoding="utf-8")
        exp_dirs.append(exp_dir)
        
    print(f"✅ {len(experiments)} configurations générées.")

    # 2. Découpage en lots (Chunks) pour contourner la limite SLURM
    CHUNK_SIZE = 5  # Tu peux l'augmenter à 8 ou 10 selon la limite de ton cluster
    chunks = [exp_dirs[i:i + CHUNK_SIZE] for i in range(0, len(exp_dirs), CHUNK_SIZE)]
    print(f"📦 Découpage en {len(chunks)} lots de {CHUNK_SIZE} expériences maximum.")

    # 3. Boucle sur les lots
    for i, chunk in enumerate(chunks):
        list_file = SLURM_DIR / f"experiments_list_chunk{i}.txt"
        list_file.write_text("\n".join(str(d) for d in chunk) + "\n", encoding="utf-8")
        
        task_script = SLURM_DIR / f"run_task_chunk{i}.sh"
        task_script.write_text(f"""#!/usr/bin/env bash
set -euo pipefail
EXP_DIR=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "{list_file}")
export OMP_NUM_THREADS=${{SLURM_CPUS_PER_TASK:-1}}
cd "$EXP_DIR"
"{EXECUTABLE}" "config.toml" > "run.log" 2>&1
""", encoding="utf-8")
        task_script.chmod(0o755)

        sbatch_script = SLURM_DIR / f"submit_array_chunk{i}.sbatch"
        sbatch_script.write_text(f"""#!/usr/bin/env bash
#SBATCH --job-name=Ferro_C{i}
#SBATCH --output={SLURM_DIR}/slurm_chunk{i}_%A_%a.out
#SBATCH --error={SLURM_DIR}/slurm_chunk{i}_%A_%a.err
#SBATCH --array=0-{len(chunk) - 1}
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00

bash {task_script}
""", encoding="utf-8")

        print(f"🚀 Soumission du lot {i+1}/{len(chunks)} ({len(chunk)} jobs). En attente de la fin du lot...")
        
        # Le paramètre --wait indique à sbatch de bloquer l'exécution de Python
        # jusqu'à ce que tous les jobs de cet array soient terminés.
        subprocess.run(["sbatch", "--wait", str(sbatch_script)], check=True)
        
    print("🎉 Toutes les expériences sont terminées !")

if __name__ == "__main__":
    if not EXECUTABLE.exists():
        print(f"❌ Erreur : L'exécutable {EXECUTABLE} est introuvable. Compile le projet C++ d'abord.")
    else:
        setup_and_submit_slurm()