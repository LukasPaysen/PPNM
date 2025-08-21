# conv_tol.gnu
set terminal pngcairo size 800,600 enhanced font 'Arial,14'
set output 'conv_tol.png'

set title "Convergence vs Tolerance"
set xlabel "Tolerance"
set ylabel "Computed E0"
set grid
set logscale x 10  # if tol spans orders of magnitude

plot "conv_tol.dat" using 1:2 with linespoints lw 2 pt 7 title "E0(tol)"
