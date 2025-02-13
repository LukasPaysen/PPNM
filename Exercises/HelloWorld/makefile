SRC = main.cs
OUT = main.exe

default: Out.txt
	cat $<
# Out.txt depends on the exe file
Out.txt : $(OUT)
	mono $(OUT) > Out.txt # run the exe, send output to Out.txt

# Compile the .cs file into a .exe
$(OUT) : $(SRC)
	mcs -target:exe -out:$(OUT) $(SRC) # compile .cs to .exe

# Clean up generated files
clean:
	rm -f Out.txt $(OUT) # remove output and exe files

