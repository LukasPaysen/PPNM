default: Out.txt
	cat $<
N=1e9
TIME = time --portability --append --output Out.txt
Out.txt : main.exe
	>$@
	$(TIME) mono main.exe -terms:$(N) -threads:1 >>$@              # Out.txt depends on hello.exe
	$(TIME) mono main.exe -terms:$(N) -threads:2 >>$@
	$(TIME) mono main.exe -terms:$(N) -threads:3 >>$@
	$(TIME) mono main.exe -terms:$(N) -threads:4 >>$@
	$(TIME) mono main.exe -terms:$(N) -threads:5 >>$@
	$(TIME) mono main.exe -terms:$(N) -threads:6 >>$@
	$(TIME) mono main.exe -terms:$(N) -threads:7 >>$@
	mono main.exe -terms:$(N) -threads:8 >>$@
main.exe : main.cs             # hello.exe depends on hello.cs
	mcs -target:exe -out:main.exe main.cs # compile hello.cs, save bytecode in hello.exe
clean:                           # a phoney target, no dependencies
	rm -f Out.txt main.exe  # remove secondary files

