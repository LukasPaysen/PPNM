# plot_wavefunc.gnu
set terminal pngcairo size 800,600 enhanced font 'Arial,14'
set output 'wavefunc.png'

set title "Hydrogen 1s Wavefunction"
set xlabel "r"
set ylabel "f(r)"
set grid

plot "wavefunc.dat" using 1:2 with lines lw 2 lc rgb "blue" title "Numerical",\
     "wavefunc.dat" using 1:3 with lines lw 2 lc rgb "red" title "Analytic r*exp(-r)"
