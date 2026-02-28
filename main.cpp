//g++ main.cpp -o main -std=c++17 -l ncurses -l tinfo -l gmpxx -l gmp -O3 -march=native -funsafe-loop-optimizations
#include <iostream>
#include <cmath>
#include <gmp.h>
#include <ncurses.h>

#include <fstream>
std::fstream file("out", std::ios::trunc | std::ios::out);

//tweak to suit your font, most work at 2:1
#define RATIO_H 2.0
#define RATIO_W 1.0

char* grid;
constexpr char density[] = {" .,%&@"};
constexpr int d_size = sizeof(density)/sizeof(density[0])-2;

//iterations, perturbation threshold & global stuff
double threshold = 1.0e-15;
constexpr int iters = 5000;
double dx;
double dy;
double width_d;
double height_d;

//ncurses variables, zoom constant
#define zoom 0.05
#define norm_lim 3
int row, col;
char ch;

//constants for coloring & antialiasing
constexpr float glow = 150;
constexpr float offset[] = {0, 0.25, -0.25};
constexpr double inv_ln = 1.0/std::log(norm_lim+1);

//perturbation stuff
double* store_Z = new double[(iters+1)*2];
mpf_t temp, R_H, R_W, mpf_zoom;

struct perturb_data{
        //i keep track of current reference in respect
        //to the array grid, not complex plane
        int Z_r;
        int Z_i;
        int ref_iters;
};

//complex plane bounds
struct bounds{
        mpf_t left;
        mpf_t right;
        mpf_t top;
        mpf_t bottom;

        mpf_t height;
        mpf_t width;
} bound;

struct mariani_data{
        mariani_data(const int h_o, const int w_o, const int h, const int w, const int m_i, const int m_r):
                h_offset(h_o), w_offset(w_o), height(h), width(w), midp_i(m_i), midp_r(m_r){}
        mariani_data(const int row, const int col){
                h_offset = 0;
                w_offset = 0;
                height = row;
                width = col;
                midp_i = row/2;
                midp_r = col/2;
        }
        int h_offset;
        int w_offset;
        int midp_i;
        int midp_r;
        int height;
        int width;
};

double calculate_pixel(const double c_r, const double c_i){
        double z_r = 0;
        double z_i = 0;
        double z2_r = 0;
        double z2_i = 0;
        //DEM
        double dz_r = 1;
        double dz_i = 0;
        double dz_r_new;
        double dz_i_new;
        double temp;
        //Periodicity
        float z_r_old = 0;
        float z_i_old = 0;
        int period = 0;
        for (int n=0; n<iters; n++){
                dz_r_new = 2*(z_r*dz_r-z_i*dz_i)+1;
                dz_i_new = 2*(z_r*dz_i+z_i*dz_r);
                dz_r = dz_r_new;
                dz_i = dz_i_new;

                z2_r = z_r * z_r;
                z2_i = z_i * z_i;
                z_i = (2.0*z_i*z_r)+c_i;
                z_r = z2_r-z2_i+c_r;

                temp = std::sqrt(z2_r+z2_i);
                if (temp > 4){
                        //file<<temp*std::log(temp)/std::sqrt(dz_r*dz_r+dz_i*dz_i)<<"\n";
                        return temp*std::log(temp)/std::sqrt(dz_r*dz_r+dz_i*dz_i);
                }
                if ((float(z_r)==z_r_old) && (float(z_i)==z_i_old)){
                        return 0;
                }
                period++;
                if (period >= 20){
                        period = 0;
                        z_r_old = z_r;
                        z_i_old = z_i;
                }
        }
        return 0;
}

int reference_orbit(const mpf_t c_r, const mpf_t c_i){
        mpf_t z_r;
        mpf_t z_i;
        mpf_t z2_r;
        mpf_t z2_i;
        mpf_t temp;
        mpf_inits(z_r, z_i, z2_r, z2_i, temp, NULL);
        for (int n=0; n<iters; n++){
                store_Z[n*2]=mpf_get_d(z_r);
                store_Z[n*2+1]=mpf_get_d(z_i);

                mpf_mul(z2_r, z_r, z_r);
                mpf_mul(z2_i, z_i, z_i);
                //imag
                mpf_mul(temp, z_i, z_r);
                mpf_mul_ui(temp, temp, 2);
                mpf_add(z_i, c_i, temp);
                //real
                mpf_sub(temp, z2_r, z2_i);
                mpf_add(z_r, c_r, temp);

                mpf_add(temp, z2_r, z2_i);
                if (mpf_cmp_ui(temp, 16) > 0){
                        mpf_clears(z_r, z_i, z2_r, z2_i, temp, NULL);
                        return n;
                }
        }
        mpf_clears(z_r, z_i, z2_r, z2_i, temp, NULL);
        return iters;
}

double pixel_orbit(const double d_r, const double d_i, const int ref_iters){
//this shit was so hard
#define Z_r store_Z[n*2]
#define Z_i store_Z[n*2+1]
#define Zz_r (Z_r+z_r)
#define Zz_i (Z_i+z_i)
        double z_r = 0;
        double z_i = 0;
        double z2_r = 0;
        double z2_i = 0;
        double temp;
        //DEM
        double dz_r = 1;
        double dz_i = 0;
        double dz_r_new;
        double dz_i_new;
        //Periodicity
        float z_r_old = 0;
        float z_i_old = 0;
        int period = 0;
        for (int n=0; n<ref_iters; n++){
                dz_r_new = 2*(Zz_r*dz_r-Zz_i*dz_i)+1;
                dz_i_new = 2*(Zz_r*dz_i+Zz_i*dz_r);
                dz_r = dz_r_new;
                dz_i = dz_i_new;

                z2_r = z_r * z_r;
                z2_i = z_i * z_i;

                temp = z_r;
                z_r = 2*(Z_r*z_r-Z_i*z_i)+z2_r-z2_i+d_r;
                z_i = 2*(Z_r*z_i+Z_i*temp+temp*z_i)+d_i;
                temp = std::sqrt(Zz_r*Zz_r+Zz_i*Zz_i);

                if (temp > 4.0){
                        dz_r_new = 2*(z_r*dz_r-z_i*dz_i)+1;
                        dz_i_new = 2*(z_r*dz_i+z_i*dz_r);
                        dz_r = dz_r_new;
                        dz_i = dz_i_new;
                        return temp*std::log(temp)/std::sqrt(dz_r*dz_r+dz_i*dz_i);
                }
                if ((float(z_r)==z_r_old) && (float(z_i)==z_i_old)){
                        return 0;
                }
                period++;
                if (period >= 20){
                        period = 0;
                        z_r_old = z_r;
                        z_i_old = z_i;
                }
        }
        return 0;
#undef Z_r
#undef Z_i
#undef Zz_r
#undef Zz_i
}

int get_ref() {
        return 0;
}

template <bool perturb>
char process_pixel(const double& dy, const double& dx, const double& c_real, const double& c_imag, const perturb_data* const p_d){
        double modulus;
        int cnt = 0;
        double mu = 0;
        //antialiasing yaay
        for (int i=0; i<2; i++){
                double offset_real = c_real+offset[i]*dy;
                for (int j=0; j<2; j++){
                        double offset_imag = c_imag+offset[j]*dx;
                        if constexpr(perturb){
                                //technically c here would be delta
                                modulus = pixel_orbit(offset_real, offset_imag, p_d->ref_iters);
                        }
                        else {
                                modulus = calculate_pixel(offset_real, offset_imag);
                        }
                        if (modulus > 0){
                                cnt++;
                                mu+=modulus;
                        }
                }
        }
        if (cnt > 1){
                double norm = (mu*glow)/(cnt*dx);
                if (norm > norm_lim) norm = norm_lim;
                double index = std::log(1.0+norm)*inv_ln;
                //this +0.1 works on hopes and dreams:
                //i hope and dream it does not segfault
                return density[int((1.0-index)*(d_size+0.1))];
                //it probably mathematically can't? idk
        }
        else {
                return density[0];
        }
}

bool check_perimeter(mariani_data m_d){
        char ch = grid[m_d.h_offset*col+m_d.w_offset];

        if (m_d.h_offset+m_d.height >= row){
                m_d.height = row-m_d.h_offset-1;
        }
        if (m_d.w_offset+m_d.width >= col){
                m_d.width = col-m_d.w_offset-1;
        }

        for (int x=0; x<m_d.height; x++){
                if (grid[(m_d.h_offset+x)*col+(m_d.w_offset+m_d.width)] != ch){
                        return 0;
                }
                if (grid[(m_d.h_offset+x)*col+m_d.w_offset] != ch){
                        return 0;
                }
        }

        for (int y=0; y<m_d.width; y++){
                if (grid[m_d.h_offset*col+(m_d.w_offset+y)] != ch){
                        return 0;
                }
                if (grid[(m_d.h_offset+m_d.height)*col+(m_d.w_offset+y)] != ch){
                        return 0;
                }
        }
        return 1;
}

template <bool perturb>
void initial_perimeter(const double& dy, const double& dx, const mariani_data& m_d, const perturb_data* const p_d){
        if constexpr(perturb){
                double d_real;
                double d_imag;
                double temp;
                for (int x=0; x<row; x++){
                        temp = x-m_d.midp_i;
                        d_imag = temp*dx;

                        temp = -m_d.midp_r;
                        d_real = temp*dy;
                        grid[x*col] = process_pixel<perturb>(dy, dx, d_real, d_imag, p_d);
                        temp = m_d.midp_r;
                        d_real = temp*dy;
                        grid[x*col+col-1] = process_pixel<perturb>(dy, dx, d_real, d_imag, p_d);
                }
                for (int y=0; y<col; y++){
                        temp = y-m_d.midp_r;
                        d_real= temp*dy;

                        temp = -m_d.midp_i;
                        d_imag = temp*dx;
                        grid[y] = process_pixel<perturb>(dy, dx, d_real, d_imag, p_d);
                        temp = m_d.midp_i;
                        d_imag = temp*dx;
                        grid[(row-1)*col+y] = process_pixel<perturb>(dy, dx, d_real, d_imag, p_d);
                }
        }
        else {

                double c_real;
                double c_imag;
                double left_d = mpf_get_d(bound.left);
                double bottom_d = mpf_get_d(bound.bottom);

                for (int x=0; x<row; x++){
                        c_imag = bottom_d+(x*dx);

                        c_real = left_d+dy;
                        grid[x*col] = process_pixel<perturb>(dy, dx, c_real, c_imag, p_d);
                        c_real = left_d+((col-1)*dy);
                        grid[x*col+col-1] = process_pixel<perturb>(dy, dx, c_real, c_imag, p_d);
                }

                for (int y=0; y<col; y++){
                        c_real = left_d+(y*dy);

                        c_imag = bottom_d+dx;
                        grid[y] = process_pixel<perturb>(dy, dx, c_real, c_imag, p_d);
                        c_imag = bottom_d+((row-1)*dx);
                        grid[(row-1)*col+y] = process_pixel<perturb>(dy, dx, c_real, c_imag, p_d);
                }
        }
}

template <bool perturb>
void mariani_silver(const double& dy, const double& dx, const mariani_data m_d, const perturb_data*  const p_d){
        const int round_h = m_d.height/2;
        const int round_w = m_d.width/2;
        const int midp_h = m_d.height-round_h;
        const int midp_w = m_d.width-round_w;

        if (check_perimeter(m_d) == 1 && row/1.5 > m_d.height) {
                for (int y=0; y<m_d.width-0; y++){
                        for (int x=0; x<m_d.height-0; x++){
                                grid[(m_d.h_offset+x)*col+(m_d.w_offset+y)] = grid[m_d.h_offset*col+m_d.w_offset];
                        }
                }
                return;
        }

        if constexpr(perturb){
                //d is for delta btw (btw)
                double d_real;
                double d_imag;
                double temp;
                if (m_d.height <= 4 || m_d.width <= 4){
                        for (int y=0; y<m_d.width; y++){
                                temp = (y+m_d.w_offset)-m_d.midp_r;
                                d_real = temp*dy;
                                for (int x=0; x<m_d.height;x++){
                                        temp = (x+m_d.h_offset)-m_d.midp_i;
                                        d_imag = temp*dx;
                                        grid[(m_d.h_offset+x)*col+(m_d.w_offset+y)] = process_pixel<perturb>(dy, dx, d_real, d_imag, p_d);
                                }
                        }
                        return;
                }

                temp = (m_d.w_offset+round_w)-m_d.midp_r;
                d_real = temp*dy;
                for (int x=0; x<m_d.height; x++){
                        temp = (x+m_d.h_offset)-m_d.midp_i;
                        d_imag = temp*dx;
                        grid[(m_d.h_offset+x)*col+(m_d.w_offset+round_w)] = process_pixel<perturb>(dy, dx, d_real, d_imag, p_d);
                }

                temp = (m_d.h_offset+round_h)-m_d.midp_i;
                d_imag = temp*dx;
                for (int y=0; y<m_d.width; y++){
                        temp = (y+m_d.w_offset)-m_d.midp_r;
                        d_real = temp*dy;
                        grid[(m_d.h_offset+round_h)*col+(m_d.w_offset+y)] = process_pixel<perturb>(dy, dx, d_real, d_imag, p_d);
                }
        }

        else {
                double c_real;
                double c_imag;
                double left_d = mpf_get_d(bound.left);
                double bottom_d = mpf_get_d(bound.bottom);

                if (m_d.height <= 4 || m_d.width <= 4){
                        for (int y=0; y<m_d.width; y++){
                                c_real = left_d+((m_d.w_offset+y)*dy);
                                for (int x=0; x<m_d.height; x++){
                                        c_imag = bottom_d+((m_d.h_offset+x)*dx);
                                        grid[(m_d.h_offset+x)*col+(m_d.w_offset+y)] = process_pixel<perturb>(dy, dx, c_real, c_imag, p_d);
                                }
                        }
                        return;
                }
                c_real = left_d+((m_d.w_offset+round_w)*dy);
                for (int x=0; x<m_d.height; x++){
                        c_imag = bottom_d+((m_d.h_offset+x)*dx);
                        grid[(m_d.h_offset+x)*col+(m_d.w_offset+round_w)] = process_pixel<perturb>(dy, dx, c_real, c_imag, p_d);
                }
                c_imag = bottom_d+((m_d.h_offset+round_h)*dx);
                for (int y=0; y<m_d.width; y++){
                        c_real = left_d+((m_d.w_offset+y)*dy);
                        grid[(m_d.h_offset+round_h)*col+(m_d.w_offset+y)] = process_pixel<perturb>(dy, dx, c_real, c_imag, p_d);
                }
        }

        //call 4 new rectangles
        mariani_silver<perturb>(dy, dx, {m_d.h_offset        , m_d.w_offset        , round_h, round_w, m_d.midp_i, m_d.midp_r}, p_d);
        mariani_silver<perturb>(dy, dx, {m_d.h_offset        , m_d.w_offset+round_w, round_h, midp_w , m_d.midp_i, m_d.midp_r}, p_d);
        mariani_silver<perturb>(dy, dx, {m_d.h_offset+round_h, m_d.w_offset        , midp_h , round_w, m_d.midp_i, m_d.midp_r}, p_d);
        mariani_silver<perturb>(dy, dx, {m_d.h_offset+round_h, m_d.w_offset+round_w, midp_h , midp_w, m_d.midp_i, m_d.midp_r }, p_d);
        return;
}


void generate_set(){
        double width_d = mpf_get_d(bound.width);
        double height_d = mpf_get_d(bound.height);
        const double dy=(width_d)/col;
        const double dx=(height_d)/row;

        mariani_data m_d(row, col);

        //check to perturb the set
        if (dx < threshold) {
                perturb_data p_d;
                mpf_t ref_r;
                mpf_t ref_i;
                mpf_inits(ref_r, ref_i, NULL);

                mpf_set_d(temp, m_d.midp_r*dy);
                mpf_add(ref_r, bound.left, temp);

                mpf_set_d(temp, m_d.midp_i*dx);
                mpf_add(ref_i, bound.bottom, temp);

                p_d.ref_iters = reference_orbit(ref_r, ref_i);
                initial_perimeter<1>(dy, dx, m_d, &p_d);
                mariani_silver<1>(dy, dx, m_d, &p_d);
        }
        else {
                initial_perimeter<0>(dy, dx, m_d, nullptr);
                mariani_silver<0>(dy, dx, m_d, nullptr);
        }
}

void output(){
        for (int x=0; x<row; x++){
                mvaddnstr(x, 0, &grid[x * col], col);
        }
        refresh();
}

void update_precision(const double& dx){
        //dont know why log2, internet says so
        int req_bits = int(std::abs(std::log2(dx))+32);
        if (req_bits > mpf_get_default_prec()){
                mpf_set_default_prec(req_bits);
                mpf_t r_t, l_t, t_t, b_t;
                mpf_inits(r_t, l_t, t_t, b_t, NULL);

                mpf_set(r_t,  bound.right);
                mpf_set(l_t,   bound.left);
                mpf_set(t_t,    bound.top);
                mpf_set(b_t, bound.bottom);

                mpf_clears(bound.right, bound.left, bound.top, bound.bottom, bound.width, bound.height, temp, NULL);
                mpf_inits(bound.right, bound.left, bound.top, bound.bottom, bound.width, bound.height, temp, NULL);

                mpf_set(bound.right,  r_t);
                mpf_set(bound.left,   l_t);
                mpf_set(bound.top,    t_t);
                mpf_set(bound.bottom, b_t);

                mpf_clears(l_t, r_t, t_t, b_t, NULL);
        }
}


int main (int argc, char* argv[]){

        //start curses
        initscr();
        raw();
        noecho();
        keypad(stdscr, TRUE);
        scrollok(stdscr, FALSE);
        nodelay(stdscr, TRUE);

        getmaxyx(stdscr, row, col);

        //GMP inits
        mpf_set_default_prec(64);
        mpf_inits(bound.right, bound.left, bound.top, bound.bottom, bound.width, bound.height, mpf_zoom, temp, R_H, R_W, NULL);
        mpf_set_d(R_H, RATIO_H);
        mpf_set_d(R_W, RATIO_W);
        mpf_set_d(mpf_zoom, zoom);

        //right
        mpf_set_d(bound.right, 1.5);
        //top
        mpf_mul(temp, bound.right, R_H);
        mpf_div(temp, temp, R_W);
        mpf_div_ui(temp, temp, col);
        mpf_mul_ui(bound.top, temp, row);
        //bottom & left
        mpf_neg(bound.bottom, bound.top);
        mpf_neg(bound.left, bound.right);

        grid = new char[row*col];
        for (int i=0; i<row*col; i++){
                grid[i] = ' ';
        }

        /*
        //temporary single point zoom module
        long double zoom_r=-1.7891690186048231066744683411888387638173618368159070155822017397181006156270275749142369245820396054406395755675312183271534128923049471434097690222315419202715383264050159131947029173677395015878767362;
        long double zoom_i=-0.0000003393685157671825660282302661468127283482188945938569013974696942388736569110136147219176174266842972236468548795141989046122450230790250469659063534132826324119846599278074403635939133245821264540;

        bound.left = zoom_i-2;
        bound.right = zoom_i+2;

        bound.top = zoom_r+3;
        bound.bottom = zoom_r-3;
        */

        mpf_sub(bound.width, bound.right, bound.left);
        mpf_sub(bound.height, bound.top, bound.bottom);
        width_d = mpf_get_d(bound.width);
        height_d = mpf_get_d(bound.height);
        dx=(height_d)/row;
        generate_set();
        output();

start:
        while ((ch = getch()) != 'q'){
                width_d = mpf_get_d(bound.width);
                height_d = mpf_get_d(bound.height);
                dx=(height_d)/row;
                if (ch=='h'){
                        mpf_mul(temp, bound.width, mpf_zoom);

                        mpf_sub(bound.left, bound.left, temp);
                        mpf_sub(bound.right, bound.right, temp);

                        generate_set();
                        output();
                }
                if (ch =='j'){
                        mpf_mul(temp, bound.height, mpf_zoom);

                        mpf_add(bound.top, bound.top, temp);
                        mpf_add(bound.bottom, bound.bottom, temp);

                        generate_set();
                        output();
                }
                if (ch=='k'){
                        mpf_mul(temp, bound.height, mpf_zoom);

                        mpf_sub(bound.top, bound.top, temp);
                        mpf_sub(bound.bottom, bound.bottom, temp);

                        generate_set();
                        output();
                }
                if (ch=='l'){
                        mpf_mul(temp, bound.width, mpf_zoom);

                        mpf_add(bound.left, bound.left, temp);
                        mpf_add(bound.right, bound.right, temp);

                        generate_set();
                        output();
                }
                if (ch=='='){
                        update_precision(dx);

                        mpf_mul(temp, bound.height, mpf_zoom);
                        mpf_sub(bound.top, bound.top, temp);
                        mpf_add(bound.bottom, bound.bottom, temp);


                        mpf_mul(temp, bound.width, mpf_zoom);
                        mpf_add(bound.left, bound.left, temp);
                        mpf_sub(bound.right, bound.right, temp);

                        mpf_sub(bound.width, bound.right, bound.left);
                        mpf_sub(bound.height, bound.top, bound.bottom);

                        generate_set();
                        output();
                }
                if (ch=='-'){
                        //mariani would fill the entire set if zoomed out too far
                        if (width_d > 5 || height_d > 5) goto start;

                        update_precision(dx);

                        mpf_mul(temp, bound.height, mpf_zoom);
                        mpf_add(bound.top, bound.top, temp);
                        mpf_sub(bound.bottom, bound.bottom, temp);

                        mpf_mul(temp, bound.width, mpf_zoom);
                        mpf_sub(bound.left, bound.left, temp);
                        mpf_add(bound.right, bound.right, temp);

                        mpf_sub(bound.width, bound.right, bound.left);
                        mpf_sub(bound.height, bound.top, bound.bottom);

                        generate_set();
                        output();
                }
                if (ch=='w'){
                        threshold = 1.0;
                        generate_set();
                        output();
                }
                if (ch=='e'){
                        threshold = 1.0e-16;
                        generate_set();
                        output();
                }
        }
        file.close();
        delete[] grid;
        mpf_clears(bound.right, bound.left, bound.top, bound.bottom, bound.width, bound.height, mpf_zoom, temp, R_H, R_W, NULL);

        refresh();
        endwin();
        return 0;}
