# plot_wavefunctions.gp

set terminal pngcairo size 800,600 enhanced font 'Verdana,10'
set output 'wavefunctions.png'

set title 'Hydrogen s-wave functions'
set xlabel 'r (Bohr)'
set ylabel 'f(r)'
set grid

plot \
  'wave_functions.dat' using 1:2 with lines lw 2 title 'num k=0', \
  ''                   using 1:5 with lines lw 2 dt 2 title 'exact k=0', \
  ''                   using 1:3 with lines lw 1 title 'num k=1', \
  ''                   using 1:4 with lines lw 1 title 'num k=2'
