# plot_errors.gp
# Plots actual and estimated errors vs N and compares to 1/sqrt(N)
set terminal pngcairo size 800,600 enhanced font 'Arial,12'
set output 'circle_errors.png'
set logscale xy
set xlabel "Number of points N"
set ylabel "Error"
set title "Monte Carlo Integration Errors (Unit Circle)"
set grid

plot "circle.dat" using 1:4 with points pt 7 ps 1.5 lc rgb "red" title "Actual Error", \
     "circle.dat" using 1:3 with linespoints lw 2 lc rgb "blue" title "Estimated Error", \
     "circle.dat" using 1:5 with lines lw 2 lc rgb "green" title "1/sqrt(N)"

# Repeat for polynomial integral
set output 'poly_errors.png'
set title "Monte Carlo Integration Errors (Polynomial Integral)"
plot "poly.dat" using 1:4 with points pt 7 ps 1.5 lc rgb "red" title "Actual Error", \
     "poly.dat" using 1:3 with linespoints lw 2 lc rgb "blue" title "Estimated Error", \
     "poly.dat" using 1:5 with lines lw 2 lc rgb "green" title "1/sqrt(N)"
