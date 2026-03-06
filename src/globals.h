////////////////////////////////
//General control stuff
////////////////////////////////
#ifndef _GLOBLS
#define _GLOBLS

#include <iostream>
#include <cmath>
#include <fstream>

#define norm_lim 3

//render variables
int iters = 5000;
int bailout = 4;

//coloring constants
constexpr char density[] = {" .,%&@"};
constexpr float glow = 150;
constexpr float offset[] = {0, 0.25, -0.25};

//mathematical constants
constexpr double inv_ln = 1.0/std::log(norm_lim+1);
constexpr double inv_ln2 = 1.0/std::log(2.0);
constexpr int d_size = sizeof(density)/sizeof(density[0])-2;

//output array
char* grid;

//perturbation threshold & zoom speed
double threshold = 1.0e-16;
float zoom_spd = 0.05;

//distance between 2 pixels
double dx;
double dy;

//ncurses variables
int row, col;
char ch;

//where Z values go?!
double* store_Z = new double[(iters+1)*2];

mpf_t temp, mpf_zoom;

enum mode{
        explore,
        zoom
};

std::fstream binfile;

void (*generate_set_p)(const double&, const double&);
void (*update_precision_p)(const double&);
void (*navigate_p)();

char* part_i;
char* part_r;

struct frame{
        int width;
        int height;
        char* data;
};

struct bounds{
//complex plane bounds
//i.e our view point
        mpf_t left;
        mpf_t right;
        mpf_t top;
        mpf_t bottom;

        mpf_t height;
        mpf_t width;
        double width_d;
        double height_d;
} bound;

struct perturb_data{
        mpf_t Z_r;
        mpf_t Z_i;
        int iZ_r;
        int iZ_i;
        int Z_iters;
};

struct mariani_data{
        mariani_data(const int h_o, const int w_o, const int h, const int w):
                h_offset(h_o), w_offset(w_o), height(h), width(w){}
        mariani_data(const int row, const int col){
                h_offset = 0;
                w_offset = 0;
                height = row;
                width = col;
        }
        int h_offset;
        int w_offset;
        int height;
        int width;
};

char sprite(){
        static int i;
        constexpr char ref[] = {"-\\|/"};
        i++;
        if (i > sizeof(ref)/sizeof(ref[0])-2) i = 0;
        return ref[i];
}
#endif
