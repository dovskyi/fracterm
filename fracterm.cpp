/*
        fracterm v0.2
        README for details
*/
#include <iostream>
#include <cmath>
#include <gmp.h>
#include <ncurses.h>
#include <string.h>

#include "globals.h"
#include "fractals.h"
#include "pipeline.h"

//tweak to suit your font, most work at 2:1
#define R_FONT 2.0

template<formula f, color c>
void generate_set(const double& dy, const double& dx){

        mariani_data m_d(row, col);

        if (dx < threshold) {
                perturb_data p_d;
                mpf_t temp;
                mpf_inits(p_d.Z_r, p_d.Z_i, temp, NULL);

                //real
                mpf_add(temp, bound.right, bound.left);
                mpf_div_ui(p_d.Z_r, temp, 2);

                //imag
                mpf_add(temp, bound.top, bound.bottom);
                mpf_div_ui(p_d.Z_i, temp, 2);

                //position in array
                mpf_sub(temp, bound.right, p_d.Z_r);
                p_d.iZ_r= mpf_get_d(temp)/dy;
                mpf_sub(temp, bound.top, p_d.Z_i);
                p_d.iZ_i = mpf_get_d(temp)/dx;

                p_d.Z_iters = reference_orbit<f>(p_d.Z_r, p_d.Z_i);
                initial_perimeter<1, f, c>(dy, dx, m_d, &p_d);
                mariani_silver<1, f, c>(dy, dx, m_d, &p_d);
                mpf_clear(temp);
                grid[p_d.iZ_i*col+p_d.iZ_r] = sprite();
        }
        else {
                initial_perimeter<0, f, c>(dy, dx, m_d, nullptr);
                mariani_silver<0, f, c>(dy, dx, m_d, nullptr);
        }
}

void output(){
        for (int x=0; x<row; x++){
                mvaddnstr(x, 0, &grid[x * col], col);
        }
        refresh();
}

void assign_bounds(const char* part_r, const char* part_i, mpf_t height, mpf_t width){
        mpf_t zoom_r, zoom_i, temp;
        mpf_inits(temp, zoom_i, zoom_r, NULL);
        mpf_set_str(zoom_r, part_r, 10);
        mpf_set_str(zoom_i, part_i, 10);
        if (col>=row){
                mpf_add(bound.top,    zoom_i, width);
                mpf_sub(bound.bottom, zoom_i, width);
                mpf_sub(temp, bound.top, bound.bottom);
                double dx = mpf_get_d(temp)/row;
                mpf_set_d(temp, col/2.0*dx/R_FONT);
                mpf_add(bound.right, zoom_r, temp);
                mpf_sub(bound.left,  zoom_r, temp);

        }
        else{
                mpf_add(bound.right, zoom_r, height);
                mpf_sub(bound.left,  zoom_r, height);
                mpf_sub(temp, bound.right, bound.left);
                double dy = mpf_get_d(temp)/col;
                mpf_set_d(temp, col/2.0*dy*R_FONT);
                mpf_add(bound.top,     zoom_i, temp);
                mpf_sub(bound.bottom,  zoom_i, temp);
        }
        mpf_clears(temp, zoom_r, zoom_i, NULL);
}

template <mode m>
void update_precision(const double& dx){
        //dont know why log2, internet says so
        int req_bits = int(std::abs(std::log2(dx))+32);
        if (req_bits > mpf_get_default_prec()){
                mpf_set_default_prec(req_bits);
                mpf_t r_t, l_t, t_t, b_t, w_t, h_t;
                mpf_inits(r_t, l_t, t_t, b_t, w_t, h_t, NULL);

                mpf_set(r_t,  bound.right);
                mpf_set(l_t,   bound.left);
                mpf_set(t_t,    bound.top);
                mpf_set(b_t, bound.bottom);
                mpf_set(w_t, bound.width);
                mpf_set(h_t, bound.height);

                if constexpr(m==zoom){
                        mpf_div_ui(w_t, w_t, 2);
                        mpf_div_ui(h_t, h_t, 2);
                }

                mpf_clears(bound.right, bound.left, bound.top, bound.bottom, bound.width, bound.height, temp, NULL);
                mpf_inits(bound.right, bound.left, bound.top, bound.bottom, bound.width, bound.height, temp, NULL);

                mpf_set(bound.right,  r_t);
                mpf_set(bound.left,   l_t);
                mpf_set(bound.top,    t_t);
                mpf_set(bound.bottom, b_t);

                if constexpr(m==zoom){
                        assign_bounds(part_r, part_i, h_t, w_t);
                }

                mpf_clears(l_t, r_t, t_t, b_t, w_t, h_t, NULL);
        }
}


void update_dydx(){
        bound.width_d = mpf_get_d(bound.width);
        bound.height_d = mpf_get_d(bound.height);
        dy=(bound.width_d)/col;
        dx=(bound.height_d)/row;
}

void usage_message(){
        std::cout<<"\nFracterm: fractal explorer for the terminal\n\n"
                <<"fracterm [flags]\n\n"
                <<"Flags:\n"
                <<"-h | display this message\n"
                <<"-d | set all defaults\n"
                <<"-f | set fractal [default:mandelbrot]\n"
                <<"    <mandelbrot, burning_ship, custom_formula>[no perturbation for ship yet]\n"
                <<"-c | set color   [default:DEM]\n"
                <<"    <DEM, dwell, custom_color>\n"
                <<"-i | set iterations\n"
                <<"-b | set bailout\n"
                <<"-m | set mode    [default:explore]\n"
                <<"    <explore>\n"
                <<"    <zoom [real] [imag]>\n\n"
                <<"While in terminal window:\n"
                <<"h/l | left/right\n"
                <<"j/k | down/up\n"
                <<"-/= | zoom in/out\n"
                <<"q   | quit\n\n"
                <<"For detailed description of flags, visit\n"
                <<"misc/flags\n";
}

void error_message(){
        std::cout<<"\nIncorrect Arguments\n\n"
                <<"-h for help\n\n";
}

int main (int argc, char* argv[]){
        void (*generate_set_p)(const double&, const double&);
        void (*update_precision_p)(const double&);

        const char *formula_choice = "mandelbrot";
        const char *color_choice = "DEM";
        const char *mode_choice = "explore";

        if (argc <= 1){
                std::cout<<"\nFracterm: fractal explorer for the terminal\n\n"
                        <<"-d for all defaults\n"
                        <<"-h for help\n";
                return 0;
        }

        for (int i=1; i<argc; i++){
                if (argv[i][0]=='-'){
                        switch(argv[i][1]){
                                case 'd':                             break;
                                case 'f': formula_choice = argv[i+1]; break;
                                case 'c': color_choice = argv[i+1];   break;
                                case 'i': iters = atoi(argv[i+1]);    break;
                                case 'b': bailout = atoi(argv[i+1]);  break;
                                case 'h': usage_message();         return 0;
                                case 'm': mode_choice = argv[i+1];
                                          part_r = new char[sizeof(argv[i+2])];
                                          part_i = new char[sizeof(argv[i+3])];
                                          part_r = argv[i+2];
                                          part_i = argv[i+3];
                                          i+=4;
                                          break;
                                default:
                                          error_message();
                                          return 0;
                        }
                }
        }//go go if else chain choice assigner!
        if (!strcmp(formula_choice, "mandelbrot")){
                if (!strcmp(color_choice, "DEM")){
                        generate_set_p = &generate_set<mandelbrot, DEM>;
                }
                else if (!strcmp(color_choice, "dwell")){
                        generate_set_p = &generate_set<mandelbrot, dwell>;
                }
                else {
                        error_message();
                        return 0;
                }
        }
        else if (!strcmp(formula_choice, "burning_ship")){
                if (!strcmp(color_choice, "DEM")){
                        generate_set_p = &generate_set<burning_ship, DEM>;
                }
                else if (!strcmp(color_choice, "dwell")){
                        generate_set_p = &generate_set<burning_ship, dwell>;
                }
                else {
                        error_message();
                        return 0;
                }
        }
        else {
                error_message();
                return 0;
        }
        if (!strcmp(mode_choice, "explore")){
                update_precision_p = &update_precision<explore>;
        }
        else if (!strcmp(mode_choice, "zoom")){
                update_precision_p = &update_precision<zoom>;
        }
        else {
                error_message();
                return 0;
        }

        //start curses
        initscr();
        raw();
        noecho();
        keypad(stdscr, TRUE);
        scrollok(stdscr, FALSE);

        getmaxyx(stdscr, row, col);

        //GMP inits
        mpf_set_default_prec(64);
        mpf_inits(bound.right, bound.left, bound.top, bound.bottom, bound.width, bound.height, mpf_zoom, temp, NULL);
        mpf_set_d(mpf_zoom, zoom_spd);

        grid = new char[row*col];
        for (int i=0; i<row*col; i++){
                grid[i] = ' ';
        }

        mpf_t two, one;
        mpf_inits(two, one, NULL);
        mpf_set_ui(two, 2);
        mpf_set_ui(one, 1);
        if (!strcmp(mode_choice, "explore")){
                assign_bounds("0","0", two, one);
        }
        else{
                assign_bounds(part_r, part_i, two, one);
        }
        mpf_clears(two, one, NULL);

        mpf_sub(bound.width, bound.right, bound.left);
        mpf_sub(bound.height, bound.top, bound.bottom);
        bound.width_d = mpf_get_d(bound.width);
        bound.height_d = mpf_get_d(bound.height);
        dy=(bound.width_d)/col;
        dx=(bound.height_d)/row;
        generate_set_p(dy, dx);
        output();

        while ((ch = getch()) != 'q'){
                if (ch=='h'){
                        mpf_mul(temp, bound.width, mpf_zoom);

                        mpf_sub(bound.left, bound.left, temp);
                        mpf_sub(bound.right, bound.right, temp);

                        generate_set_p(dy, dx);
                        output();
                }
                if (ch =='j'){
                        mpf_mul(temp, bound.height, mpf_zoom);

                        mpf_add(bound.top, bound.top, temp);
                        mpf_add(bound.bottom, bound.bottom, temp);

                        generate_set_p(dy, dx);
                        output();
                }
                if (ch=='k'){
                        mpf_mul(temp, bound.height, mpf_zoom);

                        mpf_sub(bound.top, bound.top, temp);
                        mpf_sub(bound.bottom, bound.bottom, temp);

                        generate_set_p(dy, dx);
                        output();
                }
                if (ch=='l'){
                        mpf_mul(temp, bound.width, mpf_zoom);

                        mpf_add(bound.left, bound.left, temp);
                        mpf_add(bound.right, bound.right, temp);

                        generate_set_p(dy, dx);
                        output();
                }
                if (ch=='='){
                        update_precision_p(dx);

                        mpf_mul(temp, bound.height, mpf_zoom);
                        mpf_sub(bound.top, bound.top, temp);
                        mpf_add(bound.bottom, bound.bottom, temp);


                        mpf_mul(temp, bound.width, mpf_zoom);
                        mpf_add(bound.left, bound.left, temp);
                        mpf_sub(bound.right, bound.right, temp);

                        mpf_sub(bound.width, bound.right, bound.left);
                        mpf_sub(bound.height, bound.top, bound.bottom);

                        update_dydx();
                        generate_set_p(dy, dx);
                        output();
                }
                if (ch=='-'){
                        update_precision_p(dx);

                        mpf_mul(temp, bound.height, mpf_zoom);
                        mpf_add(bound.top, bound.top, temp);
                        mpf_sub(bound.bottom, bound.bottom, temp);

                        mpf_mul(temp, bound.width, mpf_zoom);
                        mpf_sub(bound.left, bound.left, temp);
                        mpf_add(bound.right, bound.right, temp);

                        mpf_sub(bound.width, bound.right, bound.left);
                        mpf_sub(bound.height, bound.top, bound.bottom);

                        update_dydx();
                        generate_set_p(dy, dx);
                        output();
                }
        }
        delete[] grid;
        delete[] store_Z;
        mpf_clears(bound.right, bound.left, bound.top, bound.bottom, bound.width, bound.height, mpf_zoom, temp, NULL);

        refresh();
        endwin();
        return 0;}
