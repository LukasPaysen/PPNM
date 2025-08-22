# ann_all.gp
# Generates ann_fit.png and ann_residuals.png in one run.

g(x) = cos(5*x - 1.0) * exp(-x*x)

# ---- main fit ----
set term pngcairo size 900,600
set out 'ann_fit.png'
set title "Task A: ANN fit to g(x)=cos(5x-1)·exp(-x^2)"
set xlabel 'x'
set ylabel 'y'
set grid lw 1
set key left top
plot 'ann_data.dat' using 1:2 with points pt 7 ps 1.2 title 'samples', \
     'ann_fit.dat'  using 1:2 with lines  lw 2   title 'ANN', \
     g(x) with lines lw 2 dt 2 title 'true g(x)'

# ---- residuals ----
set term pngcairo size 900,400
set out 'ann_residuals.png'
set title 'Task A: Residuals (ANN - g)'
set xlabel 'x'
set ylabel 'F_p(x) - g(x)'
set grid lw 1
set key left top
plot 'ann_fit.dat' using 1:($2 - g($1)) with lines lw 2 title 'residual (dense grid)'

unset out
