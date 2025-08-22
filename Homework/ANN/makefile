.PHONY: default plots clean

# default target: run program (stream to terminal + save) then make plots
default: Out.txt plots

Out.txt: main.exe               # run and stream output while saving to file
	mono main.exe | tee $@

# plots run after the program has generated the .dat/.gp files
plots:
	- gnuplot ann_plot.gp
	- test -f ann_residuals.gp && gnuplot ann_residuals.gp || true
	- test -f ann_all.gp && gnuplot ann_all.gp || true

main.exe: main.cs
	mcs -target:exe -out:$@ $<

clean:
	rm -f Out.txt main.exe *.png *.dat
