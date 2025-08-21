# conv_rmin.gnu
set terminal pngcairo size 800,600 enhanced font 'Arial,14'
set output 'conv_rmin.png'

set title "Convergence vs rmin"
set xlabel "rmin"
set ylabel "Computed E0"
set grid

plot "conv_rmin.dat" using 1:2 with linespoints lw 2 pt 7 title "E0(rmin)"
