# plot_qspline.gp
set terminal pngcairo size 1000,700 enhanced font 'Arial,10'
set output 'qspline.png'

set multiplot layout 2,1 title 'Quadratic spline and its antiderivative'

set xlabel 'x'
set ylabel 'y'
plot 'qspline.dat' using 1:2 with lines title 'quadratic spline', \
     'qs_data_points.dat' using 1:2 with points pointtype 7 title 'data points'

set ylabel 'integral'
plot 'qspline.dat' using 1:3 with lines title 'qspline integral', \
     'qspline.dat' using 1:4 with lines title 'true integral (for sin)'

unset multiplot
set output
