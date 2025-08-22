set term pngcairo size 900,600
set out 'ann_fit.png'
set title 'Task A: ANN fit to g(x)=cos(5x-1)exp(-x^2) with Gaussian wavelets'
set xlabel 'x'
set ylabel 'y'
set key left top
plot \
  'ann_data.dat' using 1:2 with points pt 7 ps 1.2 title 'data', \
  'ann_fit.dat'  using 1:2 with lines  lw 2   title 'ANN approximation'
