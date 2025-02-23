default: Out.txt fig2.svg
	cat Out.txt
fig2.svg: plt.gpi Out.txt
	gnuplot plt.gpi
Out.txt : main.exe              # Out.txt depends on hello.exe
	mono main.exe > Out.txt # run hello.exe, send output to Out.txt
main.exe : main.cs             # hello.exe depends on hello.cs
	mcs -target:exe -out:main.exe main.cs # compile hello.cs, save bytecode in hello.exe
clean:                           # a phoney target, no dependencies
	rm -f Out.txt main.exe fig2.svg # remove secondary files
