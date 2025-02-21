inputfile=input.txt
default: Out.txt
	cat $<
Out.txt : main.exe              # Out.txt depends on hello.exe
	printf "1.1\n2.2\n3.3\n4.4\n5.5\n" > $(inputfile)
	mono main.exe -input:$(inputfile) -output:Out.txt # run hello.exe, send output to Out.txt
main.exe : main.cs             # hello.exe depends on hello.cs
	mcs -target:exe -out:main.exe main.cs # compile hello.cs, save bytecode in hello.exe
clean:                           # a phoney target, no dependencies
	rm -f Out.txt main.exe $(inputfile) # remove secondary files
