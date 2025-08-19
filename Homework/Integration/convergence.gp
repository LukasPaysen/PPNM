# erf_plot.gp
# Plot erf(1) convergence: acc (x) vs |error| (y) in log-log scale
# Usage: gnuplot erf_plot.gp
# Produces: erf_convergence.png

set terminal pngcairo size 900,600 enhanced font "Verdana,12"
set output "erf_convergence.png"

set logscale xy
set grid
set key left top

set xlabel "acc (requested absolute tolerance)"
set ylabel "absolute error |erf_num - erf_ref|"
set format x '%.0e'
set format y '%.0e'

# If your data file contains a header line starting with '#', gnuplot will ignore it.
# Data columns: 1=acc, 2=abs_error, 3=computed, 4=func_evals
set xrange [1e-10:1e-1]

plot "erf_convergence.dat" using 1:2 with linespoints lw 2 pt 7 ps 1.2 title "abs error"
