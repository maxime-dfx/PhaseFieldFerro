#!/usr/bin/env python3
"""
landau_roots.py
================
Test de Niveau 0 / Niveau 1 : calcule analytiquement les états d'équilibre
(racines de d(chi+W)/dp = 0) du potentiel de Landau-Devonshire pour un jeu
de coefficients materiau donne, et les compare a la valeur stationnaire
retournee par la simulation C++.

Usage
-----
1. Renseigne les coefficients materiau ci-dessous (copie-les depuis ton
   fichier de config Datafile -- ce sont les memes noms que dans Math.cpp).
2. Lance :  python3 landau_roots.py
3. Compare les racines affichees a la valeur (Px, Py) obtenue en simulation
   apres relaxation en Niveau 0 (E=0, eps=0) ou Niveau 1 (eps impose).

Le script utilise scipy.optimize pour trouver TOUTES les racines
stationnaires depuis un grillage de points de depart (car chi(p) est un
polynome de degre 8 -- plusieurs minima/maxima/selles coexistent en general).
"""

import numpy as np
from scipy.optimize import fsolve
from itertools import product

# ============================================================
# 1. COEFFICIENTS MATERIAU -- a copier depuis ta config Datafile
# ============================================================
alpha1    = -1.0     # <-- remplace par tes valeurs reelles
alpha11   = 1.0
alpha12   = 3.0      # alpha12 > alpha11 : penalise les etats mixtes (p1,p2 != 0 tous les deux)
alpha111  = 0.0
alpha112  = 0.0
alpha1111 = 0.0
alpha1112 = 0.0
alpha1122 = 0.0

b1 = 0.0   # electrostriction (Niveau 1 seulement)
b2 = 0.0
b3 = 0.0

# Deformation imposee (Niveau 1). Laisser a zero pour le Niveau 0.
eps11 = 0.0
eps22 = 0.0
eps12 = 0.0

# Champ electrique impose (Niveau 2, optionnel)
E1 = 0.0
E2 = 0.0

# ============================================================
# 2. DEFINITIONS -- identiques a Math::chi_energy / dchi_dp1 / dW_dp1 etc.
# ============================================================

def dchi_dp1(p1, p2):
    p1_2, p1_3, p1_4, p1_5, p1_7 = p1**2, p1**3, p1**4, p1**5, p1**7
    p2_2, p2_4, p2_6 = p2**2, p2**4, p2**6
    return (2.0*alpha1*p1 + 4.0*alpha11*p1_3 + 2.0*alpha12*p1*p2_2
            + 6.0*alpha111*p1_5 + 4.0*alpha112*p1_3*p2_2 + 2.0*alpha112*p1*p2_4
            + 8.0*alpha1111*p1_7 + 6.0*alpha1112*p1_5*p2_2 + 2.0*alpha1112*p1*p2_6
            + 4.0*alpha1122*p1_3*p2_4)

def dchi_dp2(p1, p2):
    p1_2, p1_4, p1_6 = p1**2, p1**4, p1**6
    p2_2, p2_3, p2_4, p2_5, p2_7 = p2**2, p2**3, p2**4, p2**5, p2**7
    return (2.0*alpha1*p2 + 4.0*alpha11*p2_3 + 2.0*alpha12*p1_2*p2
            + 6.0*alpha111*p2_5 + 4.0*alpha112*p1_2*p2_3 + 2.0*alpha112*p1_4*p2
            + 8.0*alpha1111*p2_7 + 6.0*alpha1112*p1_2*p2_5 + 2.0*alpha1112*p1_6*p2
            + 4.0*alpha1122*p1_4*p2_3)

def dW_dp1(p1, p2):
    return -b1*eps11*p1 - b2*eps22*p1 - b3*(eps12 + eps12)*p2

def dW_dp2(p1, p2):
    return -b1*eps22*p2 - b2*eps11*p2 - b3*(eps12 + eps12)*p1

def stationarity(p):
    """dh/dp = dchi/dp + dW/dp - E = 0  (equilibre : force totale nulle)"""
    p1, p2 = p
    return [dchi_dp1(p1, p2) + dW_dp1(p1, p2) - E1,
            dchi_dp2(p1, p2) + dW_dp2(p1, p2) - E2]

def chi_W_value(p1, p2):
    p1_2, p1_4, p1_6, p1_8 = p1**2, p1**4, p1**6, p1**8
    p2_2, p2_4, p2_6, p2_8 = p2**2, p2**4, p2**6, p2**8
    chi = (alpha1*(p1_2+p2_2) + alpha11*(p1_4+p2_4) + alpha12*(p1_2*p2_2)
           + alpha111*(p1_6+p2_6) + alpha112*(p1_2*p2_4+p2_2*p1_4)
           + alpha1111*(p1_8+p2_8) + alpha1112*(p1_6*p2_2+p2_6*p1_2)
           + alpha1122*(p1_4*p2_4))
    W = (-0.5*b1*(eps11*p1_2+eps22*p2_2) - 0.5*b2*(eps11*p2_2+eps22*p1_2)
         - b3*(eps12+eps12)*p1*p2)
    return chi + W - E1*p1 - E2*p2

# ============================================================
# 3. RECHERCHE DES RACINES SUR UNE GRILLE DE POINTS DE DEPART
# ============================================================

def find_roots(search_range=2.0, n_grid=15, tol=1e-9):
    """Balaie une grille de points de depart et deduplique les racines trouvees."""
    starts = np.linspace(-search_range, search_range, n_grid)
    roots = []
    for p10, p20 in product(starts, starts):
        sol, info, ier, msg = fsolve(stationarity, [p10, p20], full_output=True)
        if ier != 1:
            continue
        residual = np.max(np.abs(stationarity(sol)))
        if residual > 1e-8:
            continue
        # dedupliquer
        is_new = all(np.linalg.norm(sol - r) > tol for r in roots)
        if is_new:
            roots.append(sol)
    return roots

def classify_root(p1, p2, h=1e-4):
    """Classe une racine (min/max/selle) via la matrice hessienne numerique de chi+W."""
    def f(p1_, p2_):
        return chi_W_value(p1_, p2_)
    fxx = (f(p1+h,p2) - 2*f(p1,p2) + f(p1-h,p2)) / h**2
    fyy = (f(p1,p2+h) - 2*f(p1,p2) + f(p1,p2-h)) / h**2
    fxy = (f(p1+h,p2+h) - f(p1+h,p2-h) - f(p1-h,p2+h) + f(p1-h,p2-h)) / (4*h**2)
    det = fxx*fyy - fxy**2
    if det > 0 and fxx > 0:
        return "MINIMUM (etat stable)"
    elif det > 0 and fxx < 0:
        return "MAXIMUM"
    elif det < 0:
        return "SELLE (instable)"
    else:
        return "DEGENERE (a verifier)"

# ============================================================
# 4. AFFICHAGE
# ============================================================

if __name__ == "__main__":
    print("=" * 70)
    print("Recherche des racines de d(chi+W)/dp = E  (etats d'equilibre)")
    print("=" * 70)
    print(f"Coefficients : alpha1={alpha1}, alpha11={alpha11}, alpha12={alpha12}")
    print(f"               b1={b1}, b2={b2}, b3={b3}")
    print(f"Deformation imposee : eps11={eps11}, eps22={eps22}, eps12={eps12}")
    print(f"Champ impose        : E1={E1}, E2={E2}")
    print("-" * 70)

    roots = find_roots()
    if not roots:
        print("Aucune racine trouvee -- elargis search_range ou verifie les coefficients.")
    else:
        # tri par energie croissante
        roots_sorted = sorted(roots, key=lambda r: chi_W_value(*r))
        for p1, p2 in roots_sorted:
            energy = chi_W_value(p1, p2)
            kind = classify_root(p1, p2)
            print(f"  P = ({p1:+.6f}, {p2:+.6f})   h={energy:+.6e}   {kind}")

    print("-" * 70)
    print("A comparer avec (Px, Py) obtenu en simulation (Niveau 0/1) apres")
    print("relaxation complete, sur un maillage a 1 element, E=0 impose (ou")
    print("eps impose et fige), sans gradient (domaine homogene).")
    print("Le MINIMUM de plus basse energie h est l'etat spontane attendu")
    print("si tu pars d'une polarisation initiale aleatoire proche de zero.")
    print("=" * 70)