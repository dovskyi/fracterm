# fracterm
An interactive escape-time fractal explorer made for the terminal. Graphics are very simple! It does mean that you never get the resolution above the dimensions of your terminal, however this is not the point. Noone is seeking an ascii explorer for the details; the point is the ability to have a powerful tool, rendered in a simplistic, artistic, expressionist way.

Currently supporing zooms up to 10^-18

# Install & Usage
Use -03 & march native for compilation, helps a lot with performance
```
git clone https://github.com/mykyd/fracterm
cd fracterm
g++ main.cpp -o main -std=c++17 -l ncurses -l tinfo -l gmpxx -l gmp -O3 -march=native
```
No arguments yet, just execute
```
./main
```
Window controls:
```
	k/j: up/down
    h/l: left/right
    +/-: zoom in/out
```
# Note & TD
I have a lot of fun making this project, and I will continue to add more features and optimizations until I hit my math ceiling. Most math stuff in here is based on documentation from great sources/people, translated into code. This means if I can understand it, I can code it. There are a lot of features I want to add, and I will do so periodically. 

Current TD:
Arbitrary precision & Perturbation
Multi-threading
More fractals
Data display
Video export
