# plot_eps0_vs_rmax.gp
set terminal pngcairo size 800,600 enhanced font 'Verdana,10'
set output 'epsilon0_vs_rmax.png'

set title 'Ground‐state Energy ε₀ vs. rmax'
set xlabel 'rmax (Bohr)'
set ylabel 'ε₀ (Hartree)'
set grid

# plot numeric data and the exact ε₀=–0.5 line
plot \
  'conv_rmax_vs_eps0.dat' using 1:2 with lines lw 2 title 'Numeric ε₀(rmax)', \
  -0.5                  with lines lt 2 dt 2 title 'Exact ε₀ = –0.5'
