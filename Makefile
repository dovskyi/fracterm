CXX = g++

C++FLAGS =  -std=c++17 -O3 -march=native -funsafe-loop-optimizations
CFLAGS = -std=c99 -O3 -march=native

C++LIBS = -l ncurses -l tinfo -l gmpxx -l gmp
CLIBS = -l ncurses -l tinfo

all: fracterm cinematograph

fracterm: src/fracterm.cpp
	$(CXX) src/fracterm.cpp -o fracterm $(C++FLAGS) $(C++LIBS)

cinematograph: replay/cinematograph.c
	$(CXX) replay/cinematograph.c -o cinematograph $(CFLAGS) $(CLIBS)

clean:
	rm cinematograph fracterm
