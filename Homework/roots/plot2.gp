# conv_rmax.gnu
set terminal pngcairo size 800,600 enhanced font 'Arial,14'
set output 'conv_rmax.png'

set title "Convergence vs rmax"
set xlabel "rmax"
set ylabel "Computed E0"
set grid

plot "conv_rmax.dat" using 1:2 with linespoints lw 2 pt 7 title "E0(rmax)"
