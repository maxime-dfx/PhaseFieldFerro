# ==============================================================================
# plot_diagnostics.gp
# Usage : gnuplot -e "csv='energies.csv'" plot_diagnostics.gp
#         (ou juste `gnuplot plot_diagnostics.gp` si le fichier s'appelle energies.csv
#          et se trouve dans le meme dossier)
#
# Attendu : CSV avec en-tete
# LoadStep,Time,U,W,Chi,Elec,Bulk,Surface,Total,Vmin,Vmax,SolveFailed
# ==============================================================================

if (!exists("csv")) csv = "run_job702//energies.csv"

set datafile separator ","
set key autotitle columnhead   # utilise l'entete comme legende
set grid

# Palette simple, coherente sur toutes les figures
c_bulk    = "#0F6E56"   # teal
c_surface = "#D85A30"   # coral
c_total   = "#26215C"   # purple fonce
c_U       = "#378ADD"   # blue
c_W       = "#BA7517"   # amber
c_Chi     = "#993556"   # pink
c_Elec    = "#639922"   # green
c_vmin    = "#042C53"   # blue fonce
c_vmax    = "#888780"   # gris
c_fail    = "#E24B4A"   # rouge

set terminal pngcairo size 1400,1600 enhanced font "Helvetica,11"

# ------------------------------------------------------------------
# Figure combinee (4 sous-graphes empiles, partagent l'axe x = Time)
# ------------------------------------------------------------------
set output "diagnostics_overview.png"
set multiplot layout 4,1 title "Diagnostics de simulation - {/*0.9 ".csv."}"

# --- 1) Energie totale : Bulk / Surface / Total ---
set title "Energie : Bulk vs Surface vs Total"
set xlabel ""
set ylabel "Energie (normalisee)"
unset key
set key top left
plot csv using 2:7 with linespoints pt 7 ps 0.4 lc rgb c_bulk    title "Bulk", \
     csv using 2:8 with linespoints pt 7 ps 0.4 lc rgb c_surface title "Surface", \
     csv using 2:9 with lines      lw 2         lc rgb c_total   title "Total"

# --- 2) Contributions individuelles a Bulk : U, W, Chi, Elec ---
set title "Contributions au terme Bulk"
set ylabel "Energie (normalisee)"
plot csv using 2:3 with linespoints pt 7 ps 0.4 lc rgb c_U    title "U (paroi de domaine)", \
     csv using 2:4 with linespoints pt 7 ps 0.4 lc rgb c_W    title "W (electro-elastique)", \
     csv using 2:5 with linespoints pt 7 ps 0.4 lc rgb c_Chi  title "Chi (Landau-Devonshire)", \
     csv using 2:6 with linespoints pt 7 ps 0.4 lc rgb c_Elec title "Elec (electrostatique)"

# --- 3) Champ v : Vmin / Vmax ---
set title "Champ de fracture v : min / max"
set ylabel "v"
set yrange [0:1.05]
plot csv using 2:10 with linespoints pt 7 ps 0.4 lc rgb c_vmin title "Vmin", \
     csv using 2:11 with linespoints pt 7 ps 0.4 lc rgb c_vmax title "Vmax", \
     0.02 with lines dt 2 lw 1.5 lc rgb c_fail title "seuil alpha (irreversibilite)"
set yrange [*:*]

# --- 4) Echecs de resolution (0/1) ---
set title "Echec de resolution lineaire (SolveFailed)"
set xlabel "Time"
set ylabel "0 = OK / 1 = echec"
set yrange [-0.1:1.1]
plot csv using 2:12 with impulses lw 2 lc rgb c_fail title "SolveFailed"
set yrange [*:*]

unset multiplot
unset output

# ------------------------------------------------------------------
# Figure individuelle : Vmin vs LoadStep, utile pour reperer
# precisement le pas d'amorcage/propagation de la fissure
# ------------------------------------------------------------------
set terminal pngcairo size 1200,700 enhanced font "Helvetica,12"
set output "vmin_vs_loadstep.png"
unset multiplot
set title "Evolution de Vmin en fonction du pas de charge"
set xlabel "LoadStep"
set ylabel "Vmin"
set yrange [0:1.05]
set key top right
plot csv using 1:10 with linespoints pt 7 ps 0.6 lc rgb c_vmin title "Vmin", \
     0.02 with lines dt 2 lw 1.5 lc rgb c_fail title "seuil alpha"
set yrange [*:*]
unset output

# ------------------------------------------------------------------
# Figure individuelle : energie de surface (mesure indirecte de la
# longueur de fissure, cf. Fig. 6 du papier Abdollahi & Arias 2011)
# ------------------------------------------------------------------
set output "surface_energy_vs_loadstep.png"
set title "Energie de surface (mesure indirecte de la longueur de fissure)"
set xlabel "LoadStep"
set ylabel "Surface energy"
set key top left
plot csv using 1:8 with linespoints pt 7 ps 0.6 lc rgb c_surface title "Surface"
unset output

print "Figures generees : diagnostics_overview.png, vmin_vs_loadstep.png, surface_energy_vs_loadstep.png"
