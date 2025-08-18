# plot_u_neg_u.gp - plots u and u' from u_neg_u.dat
set terminal pngcairo size 800,500 enhanced font 'Arial,12'
set output 'u_neg_u.png'
set title "Solution of u'' = -u  (u and u')"
set xlabel "t"
set grid
# columns: 1=t, 2=u, 3=u'
plot 'u_neg_u.dat' using 1:2 with lines lw 2 title 'u(t)', \
     'u_neg_u.dat' using 1:3 with lines lw 2 title "u'(t)"
