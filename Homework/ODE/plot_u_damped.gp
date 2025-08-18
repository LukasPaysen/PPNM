# plot_u_damped.gp - plots u and u' from u_damped.dat
set terminal pngcairo size 800,500 enhanced font 'Arial,12'
set output 'u_damped.png'
set title "Damped oscillator: u'' = -0.1 u' - u"
set xlabel "t"
set grid
plot 'u_damped.dat' using 1:2 with lines lw 2 title 'u(t)', \
     'u_damped.dat' using 1:3 with lines lw 2 title "u'(t)"
