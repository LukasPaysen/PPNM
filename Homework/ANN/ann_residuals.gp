# ann_residuals.gp
set term pngcairo size 900,400
set out 'ann_residuals.png'

set title 'Task A: Residuals (ANN - g)'
set xlabel 'x'
set ylabel 'F_p(x) - g(x)'
set grid lw 1
set key left top

g(x) = cos(5*x - 1.0) * exp(-x*x)

plot 'ann_fit.dat' using 1:($2 - g($1)) with lines lw 2 title 'residual (dense grid)'
