default: Out.txt
	cat $<
Out.txt : main.exe              # Out.txt depends on hello.exe
	mono main.exe | tee Out.txt
	gnuplot plot_wavefunc.gp
	gnuplot plot2.gp
	gnuplot plot3.gp
	gnuplot plot4.gp
main.exe : main.cs             # hello.exe depends on hello.cs
	mcs -target:exe -out:main.exe main.cs # compile hello.cs, save bytecode in hello.exe
clean:                           # a phoney target, no dependencies
	rm -f Out.txt main.exe  # remove secondary files
	rm -f *.png
	rm -f *.dat
