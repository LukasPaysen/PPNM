# plot.gp
# declipping.dat: 1:i  2:t  3:x_true  4:y_clipped  5:x_declipped  6:clipped_flag
set datafile separator whitespace
set term pngcairo size 1400,900
set output 'declipping.png'

set multiplot layout 2,1

# --- Top: signals ---
set title 'Signals: Clipped, Declipped, True'
set xlabel 't'
set ylabel 'Amplitude'
set key top left
plot 'declipping.dat' using 2:4 with lines lw 1 dt 3 title 'Clipped', \
     ''               using 2:5 with lines lw 2      title 'Declipped', \
     ''               using 2:3 with lines lw 2      title 'True'

# --- Bottom: error (Declipped - True only) ---
set title 'Error: Declipped - True'
set xlabel 't'
set ylabel 'Error'
set key top right
plot 'declipping.dat' using 2:($5-$3) with lines lw 2 title 'Declipped - True', \
     0 with lines lc rgb 'gray' dt 3 notitle

unset multiplot
unset output
