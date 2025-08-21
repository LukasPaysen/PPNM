set terminal pngcairo size 800,600
set output "higgs.png"
set xlabel "Energy [GeV]"
set ylabel "Signal σ(E)"
set title "Higgs Boson Breit-Wigner Fit"
set key top left
plot "data.dat" using 1:2:3 with yerrorbars title "Data", \
     "fit.dat" using 1:2 with lines lw 2 title "Breit-Wigner Fit"
