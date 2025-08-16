# plot_lin_spline.gp
set terminal pngcairo size 1000,700 enhanced font 'Arial,10'
set output 'spline.png'

set multiplot layout 2,1 title 'Linear spline and its antiderivative'

set xlabel 'x'
set ylabel 'y'
plot 'spline.dat' using 1:2 with lines title 'linear interpolant', \
     'data_points.dat' using 1:2 with points pointtype 7 title 'data points'

set ylabel 'integral'
plot 'spline.dat' using 1:3 with lines title 'spline antiderivative', \
     'spline.dat' using 1:4 with lines title 'true antiderivative sin(x)'

unset multiplot
set output

