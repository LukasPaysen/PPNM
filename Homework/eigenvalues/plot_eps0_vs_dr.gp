# plot_eps0_vs_dr.gp
set terminal pngcairo size 800,600 enhanced font 'Verdana,10'
set output 'epsilon0_vs_dr.png'

set title 'Ground‐state Energy ε₀ vs. Grid Spacing Δr'
set xlabel 'Δr (Bohr)'
set ylabel 'ε₀ (Hartree)'
set grid

# plot numeric data and the exact ε₀=–0.5 line
plot \
  'conv_dr_vs_eps0.dat' using 1:2 with lines lw 2 title 'Numeric ε₀(Δr)', \
  -0.5                  with lines lt 2 dt 2 title 'Exact ε₀ = –0.5'
