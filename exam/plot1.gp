# plot_complex.gp (same layout as before)
set datafile separator whitespace
set term pngcairo size 1400,900
set output 'declipping_complex.png'

set multiplot layout 2,1

set title 'Signals: Clipped, Declipped, True (complex)'
set xlabel 't'
set ylabel 'Amplitude'
set key top right
plot 'declipping_complex.dat' using 2:4 with lines lw 1 dt 3 title 'Clipped', \
     ''                      using 2:5 with lines lw 2      title 'Declipped', \
     ''                      using 2:3 with lines lw 2      title 'True'

set title 'Error: Declipped - True (complex)'
set xlabel 't'
set ylabel 'Error'
set key top right
plot 'declipping_complex.dat' using 2:($5-$3) with lines lw 2 title 'Declipped - True', \
     0 with lines lc rgb 'gray' dt 3 notitle

unset multiplot
unset output
