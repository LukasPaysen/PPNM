default: Out.txt

Out.txt : main.exe higgs.data.txt
	mono main.exe < higgs.data.txt | tee Out.txt
	gnuplot plot.gp
main.exe : main.cs
	mcs -target:exe -out:main.exe main.cs

clean:
	rm -f Out.txt main.exe
	rm -f *.png
	rm -f *.dat
