default: Out.txt
	cat $<
Out.txt : main.exe              # Out.txt depends on hello.exe
	echo 1 2 0 4 5 > input.txt
	mono main.exe < input.txt 2> Out.txt # run hello.exe, send output to Out.txt
main.exe : main.cs             # hello.exe depends on hello.cs
	mcs -target:exe -out:main.exe main.cs # compile hello.cs, save bytecode in hello.exe
clean:                           # a phoney target, no dependencies
	rm -f Out.txt main.exe  # remove secondary files
