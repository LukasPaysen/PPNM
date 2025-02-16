default: Out.txt
	cat $<
Out.txt : main.exe              # Out.txt depends on hello.exe
	mono main.exe > Out.txt # run hello.exe, send output to Out.txt
main.exe : main.cs sfuns.dll             # hello.exe depends on hello.cs
	mcs -r:sfuns.dll -out:main.exe main.cs # compile hello.cs, save bytecode in hello.exe
sfuns.dll: sfuns.cs  # sfuns.dll depends on sfuns.cs
	mcs -target:library -out:sfuns.dll sfuns.cs
clean:                           # a phoney target, no dependencies
	rm -f Out.txt main.exe  # remove secondary files
