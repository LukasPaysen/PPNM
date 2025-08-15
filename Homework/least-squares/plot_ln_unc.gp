set terminal pngcairo size 800,600 enhanced font 'Arial,12'
set output 'thx_ln_unc_plot.png'
set xlabel 't (days)'
set ylabel 'ln(y)'
set grid
set key left top

# Columns:
# thx_ln_data.dat:       t  ln(y)  delta_ln(y)
# thx_ln_perturbed.dat:  t  ln_central  ln_pp  ln_pm  ln_mp  ln_mm

plot 'thx_ln_data.dat' using 1:2:3 with yerrorbars pt 7 ps 1.2 title 'data ln(y)', \
     'thx_ln_perturbed.dat' using 1:2 with lines lw 2 title 'fit (central)', \
     'thx_ln_perturbed.dat' using 1:3 with lines lt 2 title '+/+', \
     'thx_ln_perturbed.dat' using 1:4 with lines lt 3 title '+/-', \
     'thx_ln_perturbed.dat' using 1:5 with lines lt 4 title '-/+', \
     'thx_ln_perturbed.dat' using 1:6 with lines lt 5 title '-/-'
