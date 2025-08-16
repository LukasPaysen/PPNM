default: Out.txt
	cat $<
Out.txt : main.exe              # Out.txt depends on hello.exe
	mono main.exe > Out.txt # run hello.exe, send output to Out.txt
	gnuplot plot_lin_spline.gp
	gnuplot plot_qspline.gp
	gnuplot plot_cubic.gp
main.exe : main.cs             # hello.exe depends on hello.cs
	mcs -target:exe -out:main.exe main.cs # compile hello.cs, save bytecode in hello.exe
clean:                           # a phoney target, no dependencies
	rm -f Out.txt main.exe  # remove secondary files
	rm -f *.png
	rm -f *.dat
