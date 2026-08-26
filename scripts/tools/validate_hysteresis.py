import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
import glob
import os

def plot_hysteresis(vtk_dir, output_image="hysteresis.png"):
    # Récupère tous les fichiers VTK générés dans l'ordre
    vtk_files = sorted(glob.glob(os.path.join(vtk_dir, "multiphysics_results*.vtk")))
    
    if not vtk_files:
        print(f"Aucun fichier VTK trouvé dans {vtk_dir}")
        return

    E_avg = []
    P_avg = []

    print(f"Analyse de {len(vtk_files)} fichiers VTK en cours...")
    for f in vtk_files:
        mesh = pv.read(f)
        
        # Les vecteurs sont exportés sous les noms "P", "E" par ton C++
        P_vectors = mesh.point_data.get("P")
        E_vectors = mesh.point_data.get("E")
        
        if P_vectors is not None and E_vectors is not None:
            # On s'intéresse à la composante Y pour l'hystérésis
            P_avg.append(np.mean(P_vectors[:, 1]))
            E_avg.append(np.mean(E_vectors[:, 1]))
            
    # Trace la boucle P-E (Polarisation vs Champ Electrique)
    plt.figure(figsize=(8, 6))
    plt.plot(E_avg, P_avg, '-o', color='#2E6FDE', linewidth=2, markersize=4)
    plt.axhline(0, color='black', linewidth=1)
    plt.axvline(0, color='black', linewidth=1)
    plt.xlabel(r"Champ Électrique Macroscopique $\bar{E}_y$ (V/m)", fontsize=12)
    plt.ylabel(r"Polarisation Macroscopique $\bar{P}_y$ (C/m$^2$)", fontsize=12)
    plt.title("Validation Expérimentale : Boucle d'Hystérésis", fontsize=14)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    
    # Sauvegarde l'image au lieu de l'afficher (parfait pour un cluster)
    plt.savefig(output_image, dpi=300)
    print(f"Graphique sauvegardé sous : {output_image}")

if __name__ == "__main__":
    dossier_vtk = "output/run_job903/VTK" 
    
    plot_hysteresis(dossier_vtk, "boucle_hysteresis.png")