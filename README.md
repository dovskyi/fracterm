# fracterm
An interactive escape-time fractal explorer made for the terminal. Graphics are very simple! It does mean that you never get the resolution above the dimensions of your terminal, however this is not the point.
Noone is seeking an ascii explorer for the details; the point is the ability to have a powerful tool, rendered in a simplistic, artistic, expressionist way.

Supporting zooms up to 1.0^-150 *[given a reference point]*

![sample](misc/pictures/sample.png)
*more images in misc/pictures*

# INSTALL
**Only supported on GNU/Linux**


Dependencies: [*dev packages*]
```
GNU Multiple Precision Arithmetic Library (GMP)
Ncurses
```
Compilation:
```
git clone https://github.com/dovskyi/fracterm
cd fracterm
make all
```
For permanent install, put the binaries [or whole repo] in */usr/local/bin*

To uninstall binaries, do *make clean*

# USAGE

```
fracterm [flags]
```
Flags:
```
-h | display this message
-d | set all defaults
-f | set fractal [default:mandelbrot]
    <mandelbrot, burning_ship, custom_formula>
    [perturbation only for mandelbrot currently]
-c | set color   [default:DEM]
    <DEM, dwell, custom_color>
-i | set iterations
-b | set bailout
-w | write all frames into a binary file
    <file name/path>
    [view using cinematograph]
-m | set mode    [default:explore]
    <explore>
    <zoom [real] [imag]>

For detailed flag description, visit
misc/docs
```
Terminal controls:
```
    k/j: up/down
    h/l: left/right
    +/-: zoom in/out
    q: quit
```
Cinematograph:
```
    cinematograph [file path] <fps>

    <fps> | An integer in range 4-48
                       [24: default]
    In window:
    a/d: forward/backward
    q  : quit
```
# Note & TD
I have a lot of fun making this project, and I will continue to add more features and optimizations until I hit my math ceiling [or get bored]. Most math stuff in here is based on documentation from great sources/people, translated into code. This means if I can understand it, I can code it. There are a lot of features I want to add, and I will do so periodically.

v0.21: I just thought it would be nice to have some sort of basic replay functionality. Now if you choose -w flag, fracterm will write every frame into a binary file. You can then view it using cinematograph [cool name, right?] written in C. /misc/docs to see how it works.

v0.2:  I finally added perturbation and DEM coloring method. Took 2 entire months to program this, but I am proud of the result. There are also a ton of new, small features; the program is unrecognizable from v0.1. The main issue right now is performance. I am hoping multithreading will fix it, we will see. I have never actually programmed multithreaded programs, so I need to do research first. But for now, shallow zooms (>1.0^-16) provide exceptional images, better than anything I could expect from ASCII graphics.

Current TD:
* Multi-threading
* Taylor series approximation
* More fractals
* Ffmpeg video export
* Ministatus, menus, etc.
