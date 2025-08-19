set terminal pngcairo size 900,650 enhanced font 'Arial,12'
set output 'strata_circle.png'
set logscale xy
set xlabel "Number of samples N"
set ylabel "Error"
set title "Stratified sampling vs Plain MC — unit circle (log-log)"
set grid
set key left top

# strata_compare_circle.dat:
# 1:N 2:plain_est 3:plain_acterr 4:strata_est 5:strata_estErr 6:strata_acterr 7:1/sqrt(N)

plot \
 'strata_compare_circle.dat' using 1:3 with points pt 7 ps 1.3 title "Plain actual error", \
 'strata_compare_circle.dat' using 1:6 with points pt 6 ps 1.3 title "Strata actual error", \
 'strata_compare_circle.dat' using 1:5 with linespoints lt 2 lw 2 pt 4 ps 1.0 title "Strata estimated error", \
 'strata_compare_circle.dat' using 1:7 with lines lt 3 lw 2 title "1/sqrt(N)"
