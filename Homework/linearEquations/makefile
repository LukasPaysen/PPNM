default: Out.txt times.data N-plot.svg

Out.txt : main.exe              # Out.txt depends on hello.exe
	mono main.exe > Out.txt # run hello.exe, send output to Out.txt
	cat $@
main.exe : main.cs             # hello.exe depends on hello.cs
	mcs -target:exe -out:main.exe main.cs # compile hello.cs, save bytecode in hello.exe
times.data: main.exe
	>$@
	for N in $$(seq 100 20 700); do \
		time --format "$$N %e" --output $@ --append mono main.exe -time -size:$$N;\
	done

N-plot.svg: times.data
	gnuplot plot.gpi



clean:                           # a phoney target, no dependencies
	rm -f Out.txt main.exe times.data N-plot.svg # remove secondary files
.PHONY: clean default
