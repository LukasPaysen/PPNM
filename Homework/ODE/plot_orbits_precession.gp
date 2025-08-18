# plot_orbits_precession.gp
set terminal pngcairo size 1000,700 enhanced font 'Arial,12'
set output 'orbits_precession.png'
set title "Orbits: circle, Newtonian ellipse, relativistic precession"
set xlabel "x"
set ylabel "y"
set size ratio -1   # equal scales in x and y
set grid

# Each data file has columns: 1=phi, 2=u(phi). Convert to (x,y) = (1/u)*cos(phi), (1/u)*sin(phi)
plot 'orbit_circle.dat'  using (1/$2)*cos($1):(1/$2)*sin($1) with lines lw 2 title "circle (eps=0, u'=0)", \
     'orbit_ellipse.dat' using (1/$2)*cos($1):(1/$2)*sin($1) with lines lw 2 title "ellipse (eps=0, u'=-0.5)", \
     'orbit_precess.dat' using (1/$2)*cos($1):(1/$2)*sin($1) with lines lw 2 title "precession (eps=0.01)"
