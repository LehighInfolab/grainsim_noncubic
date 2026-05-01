grainsim:
	g++ -O3 src/main.cpp -o grainsim.out -I$(HOME)/local/include -L$(HOME)/local/lib64 -lzip -Wl,-rpath,$(HOME)/local/lib64

# grainsim:
# 	g++ -O3 src/main.cpp -o grainsim.out
