set terminal pngcairo size 800,600 enhanced font 'Arial,12'
set output 'thx_ln_plot.png'
set xlabel 't (days)'
set ylabel 'ln(y)'
set grid
set key left top

# data file columns: t ln(y) delta_ln(y)
plot 'thx_ln_data.dat' using 1:2:3 with yerrorbars pt 7 ps 1.2 title 'data ln(y)', \
     'thx_ln_fit.dat' using 1:2 with lines lw 2 title 'fit'
