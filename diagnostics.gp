# ==============================================================================
# Script Gnuplot - PhaseFieldFerro - Comparaison avec Abdollahi & Arias (2011)
# Usage : gnuplot plot_results.gnu
# ==============================================================================

set datafile separator ","
CSV = "../output/energies.csv"   # <-- change le nom/chemin si besoin

set terminal pngcairo size 1000,750 enhanced font "Arial,12"
set grid
set key top left box

# ------------------------------------------------------------------
# 1 Energie de surface vs pas de charge (comparable a la Fig. 6)
# ------------------------------------------------------------------
set output "../results/graphiques/01_surface_energy.png"
set title "Evolution de l'energie de surface normalisee"
set xlabel "Load step"
set ylabel "Surface energy"
plot CSV using 1:8 with linespoints lw 2 pt 7 ps 0.5 lc rgb "blue" title "surface\\_energy"

# ------------------------------------------------------------------
# 2 Energie totale vs pas de charge
# ------------------------------------------------------------------
set output "../results/graphiques/02_total_energy.png"
set title "Energie totale du systeme"
set xlabel "Load step"
set ylabel "Total energy"
plot CSV using 1:9 with linespoints lw 2 pt 7 ps 0.5 lc rgb "red" title "total\\_energy"

# ------------------------------------------------------------------
# 3 Decomposition des energies (U, W, chi, electric, bulk enthalpy)
# ------------------------------------------------------------------
set output "../results/graphiques/03_energy_decomposition.png"
set title "Decomposition des contributions energetiques"
set xlabel "Load step"
set ylabel "Energy"
set key top left
plot CSV using 1:3  with lines lw 2 title "U\\_total (domain wall)", \
     CSV using 1:4  with lines lw 2 title "W\\_total (electroelastic)", \
     CSV using 1:5  with lines lw 2 title "chi\\_total (Landau-Devonshire)", \
     CSV using 1:6  with lines lw 2 title "electric\\_total", \
     CSV using 1:7  with lines lw 2 title "bulk\\_enthalpy", \
     CSV using 1:8  with lines lw 2 dt 2 title "surface\\_energy" \


# ------------------------------------------------------------------
# 4 v_min vs pas de charge : detection de l'initiation/propagation
#    (v_min < 0.02 => zone consideree comme fissuree)
# ------------------------------------------------------------------
set output "../results/graphiques/04_vmin.png"
set title "Valeur minimale du champ de phase v (fissure)"
set xlabel "Load step"
set ylabel "v\\_min"
set yrange [0:1.05]
set arrow from graph 0, first 0.02 to graph 1, first 0.02 nohead lc rgb "red" dt 2
set label "seuil de rupture v=0.02" at graph 0.02, first 0.06 tc rgb "red"
plot CSV using 1:10 with linespoints lw 2 pt 7 ps 0.5 lc rgb "black" title "v\\_min"
unset arrow
unset label

# ------------------------------------------------------------------
# 5 Longueur de fissure (proxy) vs pas de charge
# ------------------------------------------------------------------
set output "../results/graphiques/05_crack_length.png"
set title "Longueur de fissure (proxy)"
set xlabel "Load step"
set ylabel "crack\\_length\\_proxy"
set autoscale y
plot CSV using 1:12 with linespoints lw 2 pt 7 ps 0.5 lc rgb "dark-green" title "crack\\_length\\_proxy"

# ------------------------------------------------------------------
# 6 Bornes de polarisation Px, Py vs pas de charge (suivi du twinning)
# ------------------------------------------------------------------
set output "../results/graphiques/06_polarization_bounds.png"
set title "Bornes de la polarisation (suivi du twinning / switching)"
set xlabel "Load step"
set ylabel "Polarization"
set yrange [-1.2:1.2]
plot CSV using 1:14 with lines lw 2 lc rgb "blue"  title "Px\\_min", \
     CSV using 1:15 with lines lw 2 lc rgb "cyan"  title "Px\\_max", \
     CSV using 1:16 with lines lw 2 lc rgb "red"   title "Py\\_min", \
     CSV using 1:17 with lines lw 2 lc rgb "orange" title "Py\\_max"

# ------------------------------------------------------------------
# 7 Deplacements max (verif du chargement mecanique applique)
# ------------------------------------------------------------------
set output "../results/graphiques/07_displacements.png"
set title "Deplacements maximaux (verification du chargement)"
set xlabel "Load step"
set ylabel "Displacement"
unset yrange
plot CSV using 1:18 with linespoints lw 2 pt 5 ps 0.5 lc rgb "purple" title "ux\\_max\\_abs", \
     CSV using 1:19 with linespoints lw 2 pt 7 ps 0.5 lc rgb "green"  title "uy\\_max\\_abs"

print "Plots generes : 01_surface_energy.png ... 07_displacements.png"
