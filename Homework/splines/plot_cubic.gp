# plot_cubic.gp
set terminal pngcairo size 1000,700 enhanced font 'Arial,10'
set output 'cubic1.png'

set multiplot layout 2,1 title 'Cubic spline: our implementation vs gnuplot built-in'

set xlabel 'x'
set ylabel 'y'
plot 'cubic.dat' using 1:2 with lines title 'our cubic spline', \
     'cubic_nodes.dat' using 1:2 with points pointtype 7 title 'data points', \
     'cubic_nodes.dat' using 1:2 smooth csplines title 'gnuplot cspline'

set ylabel 'integral'
plot 'cubic.dat' using 1:3 with lines title 'cubic integral', \
     'cubic.dat' using 1:4 with lines title 'true integral (for sin)'

unset multiplot
set output
