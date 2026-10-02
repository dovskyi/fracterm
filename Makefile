CXX = g++

C++FLAGS =  -std=c++17 -O3 -march=native -funsafe-loop-optimizations
CFLAGS = -std=c99 -O3 -march=native

C++LIBS = $$(pkg-config --libs --cflags notcurses) $$(pkg-config --libs --cflags gmp) 

all: fracterm cinematograph

fracterm: src/fracterm.cpp
	$(CXX) src/fracterm.cpp -o fracterm $(C++FLAGS) $(C++LIBS)

clean:
	rm fracterm
