# plot_singular_compare_rand.gp
# Compare plain MC vs randomized-shift Halton QMC (mean +/- stddev)
set terminal pngcairo size 1000,700 enhanced font 'Arial,12'
set output 'singular_compare_rand.png'

set logscale xy
set xlabel "Number of samples N"
set ylabel "Error"
set title "Singular integral: plain MC vs randomized-shift Halton QMC (log-log)"
set grid
set key left top

# columns in singular_compare_rand.dat:
# 1: N
# 2: plain actual error
# 3: quasi actual error (mean across shifts)
# 4: plain estimated error (sigma*V/sqrt(N))
# 5: quasi stddev (std across shifts)  <-- use for errorbars
# 6: 1/sqrt(N)

plot \
 'singular_compare_rand.dat' using 1:2 with points pt 7 ps 1.5 title "Plain actual error", \
 'singular_compare_rand.dat' using 1:3:5 with yerrorbars pt 6 ps 1.2 title "Quasi actual error ± stddev", \
 'singular_compare_rand.dat' using 1:4 with linespoints lt 1 lw 2 pt 5 ps 1.2 title "Plain estimated error", \
 'singular_compare_rand.dat' using 1:5 with linespoints lt 2 lw 1 pt 4 ps 1.0 title "Quasi stddev (for ref)", \
 'singular_compare_rand.dat' using 1:6 with lines lt 3 lw 2 title "1/sqrt(N) (reference)"
