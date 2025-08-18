# plot_pendulum.gp - plots theta and omega from pendulum.dat
set terminal pngcairo size 800,500 enhanced font 'Arial,12'
set output 'pendulum.png'
set title "Pendulum with friction: theta(t) and omega(t)"
set xlabel "t"
set grid
plot 'pendulum.dat' using 1:2 with lines lw 2 title 'theta(t)', \
     'pendulum.dat' using 1:3 with lines lw 2 title 'omega(t)'
