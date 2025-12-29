#include <iostream>
#include <cmath>
#include <ncurses.h>

//output array & ascii color(?) palette
char* grid;
constexpr char density[] = {"    ..:#%@"};

//dynamic iterations, hard limit, & perturbation threshold
constexpr double threshold = 1.0e-40;
constexpr int iteration_lim = 10000;
int iterations = 100;

//ncurses variables, zoom constant
int row, col;
char ch;
constexpr double zoom = 0.05;

//constants for coloring & antialiasing
constexpr int K = 90;
constexpr double scale = 0.1;
constexpr double offset[] = {0, 0.25, -0.25};
constexpr long double inv_log2 = 1.0L/std::log(2.0L);

//perturbation stuff, not implemented yet
long double* store_reference = new long double[iteration_lim*2+4];
struct p_param{
        long double c_real;
        long double c_imag;
        int c_iters;
};

//"window" of the fractal
struct bounds{
        long double left = -2.0;
        long double right = 2.0;
        long double top;
        long double bottom;
};

bounds bound;

template <bool perturb>
int calculate_pixel(const long double c_r, const long double c_i, long double& modulus){
        long double z_r = 0;
        long double z_i = 0;
        long double z2_r = 0;
        long double z2_i = 0;
        //periodicity checking variables
        long double z_r_old = 0;
        long double z_i_old = 0;
        int period = 0;
        for (int n=0; n<iterations; n++){
                z2_r = z_r * z_r;
                z2_i = z_i * z_i;
                z_i = (2.0L*z_i*z_r)+c_i;
                z_r = z2_r-z2_i+c_r;
                if constexpr(perturb){
                        //z values are needed for perturbation
                        store_reference[n+1]=z_r;
                        store_reference[n+2]=z_i;
                }
                if (z2_r+z2_i > 9){
                        //extra iteration decreases the error term
                        z2_r = z_r * z_r;
                        z2_i = z_i * z_i;
                        z_i = (2.0L*z_i*z_r)+c_i;
                        z_r = z2_r-z2_i+c_r;
                        n++;
                        if constexpr(perturb){
                                store_reference[n+1]=z_r;
                                store_reference[n+2]=z_i;
                        }
                        modulus = std::sqrt(z2_r+z2_i);
                        return n;}
                if ((z_r == z_r_old) && (z_i == z_i_old)){
                        return 0;}
                period++;
                if (period >= 20){
                        period = 0;
                        z_r_old = z_r;
                        z_i_old = z_i;}
        }
        return 0;
}

/*
   int pixel_orbit(char* store_reference, const long double c0_r, const long double c0_i, const int c_iters){
   long double z0_r = 0;
   long double z0_i = 0;
   for (int n=0; n<c_iters; n++){
   z0 = 2.0L*store_reference[n]*z0 + (z0 * z0) + c0;
   if (std::norm(store_reference[n+1]+z0)>4){
   return n;}
   }
   return 0;
   }
   */

template <bool perturb>
char process_pixel(const long double& dy, const long double& dx, const long double& c0_real, const long double& c0_imag, const p_param* const p_p){

        long double modulus;
        long double mu = 0;
        int c0_iters;
        int counter = 0;
        /*Antialiasing, causes too much lag for now
        for (int i=0; i<2; i++){
                long double offsetc0_real = c0_real+offset[i]*dy;
                for (int j=0; j<2; j++){
                        long double offsetc0_imag = c0_imag+offset[j]*dx;

                        if constexpr(perturb){
                                //c0_iters = pixel_orbit(store_reference, offsetc0_real-p_param.c_real, offsetc0_imag-p_param.c_imag, p_param.c_iters);
                                //if (c0_iters > 0){
                                //      counter++;
                                //      mu = mu + (c0_iters+1 - log(log(modulus))*inv_log2);
                                //}
                        }
                        else {
                                c0_iters = calculate_pixel<0>(offsetc0_real, offsetc0_imag, modulus);
                                if (c0_iters > 0){
                                        counter++;
                                        mu = mu + (c0_iters+1 - log(log(modulus))*inv_log2);
                                }
                        }
                }
        }
        */
        c0_iters = calculate_pixel<0>(c0_real, c0_imag, modulus);
        if (c0_iters > 0){
                counter++;
                mu = mu + (c0_iters+1 - log(log(modulus))*inv_log2);
                return density[int(mu/counter*scale)%9];

        }
        /*
        if (mu > 0) {
                return density[int(mu/counter*scale)%9];
        }
        */
        else {
                return density[9];
        }
}

bool check_perimeter(const int h_offset, const int w_offset, int height, int width){
        char ch = grid[h_offset*col+w_offset];

        if (h_offset+height >= row){
                height = row-h_offset-1;
        }
        if (w_offset+width >= col){
                width = col-w_offset-1;
        }

        for (int x=0; x<height; x++){
                if (grid[(h_offset+x)*col+(w_offset+width)] != ch){
                        return 0;
                }
                if (grid[(h_offset+x)*col+w_offset] != ch){
                        return 0;
                }
        }

        for (int y=0; y<width; y++){
                if (grid[h_offset*col+(w_offset+y)] != ch){
                        return 0;
                }
                if (grid[(h_offset+height)*col+(w_offset+y)] != ch){
                        return 0;
                }
        }
        return 1;
}


template <bool perturb, bool initial_perimeter>
void mariani_silver(const long double& dy, const long double& dx, const int h_offset, const int w_offset, int height, int width, const p_param*  const p_p){
        long double c0_real;
        long double c0_imag;
        const int round_h = height/2;
        const int round_w = width/2;
        const int midp_h = height-round_h;
        const int midp_w = width-round_w;

        if constexpr(initial_perimeter){
                for (int y=0; y<width; y+=width-1){
                        c0_real = bound.bottom+(y*dy);
                        for (int x=0; x<height; x++){
                                c0_imag = bound.left+(x*dx);
                                grid[(h_offset+x)*col+(w_offset+y)] = process_pixel<perturb>(dy, dx, c0_real, c0_imag, p_p);
                        }
                }
                for (int x=0; x<height; x+=height-1){
                        c0_imag = bound.left+(x*dx);
                        for (int y=0; y<width; y++){
                                c0_real = bound.bottom+(y*dy);
                                grid[(h_offset+x)*col+(w_offset+y)] = process_pixel<perturb>(dy, dx, c0_real, c0_imag, p_p);
                        }
                }
        }

        //if perimeter has similar iterations, fill it
        if (check_perimeter(h_offset, w_offset, height, width) == 1 && row/1.5 > height) {
                for (int y=0; y<width-0; y++){
                        c0_real = bound.bottom+((w_offset+y)*dy);
                        for (int x=0; x<height-0; x++){
                                c0_imag = bound.left+((h_offset+x)*dx);
                                grid[(h_offset+x)*col+(w_offset+y)] = grid[h_offset*col+w_offset];
                        }
                }
                return;
        }

        if (height <= 4 || width <= 4){
                for (int y=0; y<width-0; y++){
                        c0_real = bound.bottom+((w_offset+y)*dy);
                        for (int x=0; x<height-0; x++){
                                c0_imag = bound.left+((h_offset+x)*dx);
                                grid[(h_offset+x)*col+(w_offset+y)] = process_pixel<perturb>(dy, dx, c0_real, c0_imag, p_p);
                        }
                }
                return;
        }
        //draw a cross to subdivide the rectangle
        c0_real = bound.bottom+((w_offset+round_w)*dy);
        for (int x=0; x<height; x++){
                c0_imag = bound.left+((h_offset+x)*dx);
                grid[(h_offset+x)*col+(w_offset+round_w)] = process_pixel<perturb>(dy, dx, c0_real, c0_imag, p_p);
        }
        c0_imag = bound.left+((h_offset+round_h)*dx);
        for (int y=0; y<width; y++){
                c0_real = bound.bottom+((w_offset+y)*dy);
                grid[(h_offset+round_h)*col+(w_offset+y)] = process_pixel<perturb>(dy, dx, c0_real, c0_imag, p_p);
        }

        //call 4 rectangles
        mariani_silver<perturb, 0>(dy, dx, h_offset        , w_offset        , round_h, round_w, p_p);
        mariani_silver<perturb, 0>(dy, dx, h_offset        , w_offset+round_w, round_h, midp_w , p_p);
        mariani_silver<perturb, 0>(dy, dx, h_offset+round_h, w_offset        , midp_h , round_w, p_p);
        mariani_silver<perturb, 0>(dy, dx, h_offset+round_h, w_offset+round_w, midp_h , midp_w , p_p);
        return;
}


void generate_set(){
        const long double dy=(bound.top-bound.bottom)/col;
        const long double dx=(bound.right-bound.left)/row;

        //check if perturbation needs to be used
        long double difference = bound.right - bound.left;

        if ((difference * difference) <= threshold) {
                p_param p_p;
                //perturbation requires a reference with long iterations
                //p_p.c_iters = best_reference(bound, p_p.c_real, p_p.c_imag);
                //if (p_p.c_iters == 0) { p_p.c_iters = iterations;}
                mariani_silver<1,1>(dy, dx, 0, 0, row, col, &p_p);
        }
        else {
                mariani_silver<0,1>(dy, dx, 0, 0, row, col, nullptr);
        }
}

void output(){
        char buffer[col+1];
        for (int x=0; x<row; x++){
                for (int y=0; y<col; y++){
                        buffer[y]=grid[x*col+y];
                }
                if (x == row-1){
                        //to not make a new line on last row
                        buffer[col-1] = '\0';
                }
                mvaddstr(x,0,buffer);
        }
        refresh();
}

int main (int argc, char* argv[]){

        //start curses
        initscr();
        raw();
        noecho();
        keypad(stdscr, TRUE);
        scrollok(stdscr, FALSE);

        getmaxyx(stdscr, row, col);
        //bounds bound(row, col);
        bound.top = (4.0L*col)/(3.0L*row);
        bound.bottom = -(4.0L*col)/(3.0L*row);

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

        generate_set();
        output();

        while ((ch = getch()) != 'q'){
                if (ch=='l'){
                        const long double width = bound.right-bound.left;
                        const long double height = bound.top-bound.bottom;
                        bound.top+=height*zoom;
                        bound.bottom+=height*zoom;
                        generate_set();
                        output();
                }
                if (ch =='h'){
                        const long double width = bound.right-bound.left;
                        const long double height = bound.top-bound.bottom;
                        bound.top-=height*zoom;
                        bound.bottom-=height*zoom;
                        generate_set();
                        output();
                }
                if (ch=='k'){
                        const long double width = bound.right-bound.left;
                        const long double height = bound.top-bound.bottom;
                        bound.left-=width*zoom;
                        bound.right-=width*zoom;
                        generate_set();
                        output();
                }
                if (ch=='j'){
                        const long double width = bound.right-bound.left;
                        const long double height = bound.top-bound.bottom;
                        bound.left+=width*zoom;
                        bound.right+=width*zoom;
                        generate_set();
                        output();
                }
                if (ch=='='){
                        const long double width = bound.right-bound.left;
                        const long double height = bound.top-bound.bottom;

                        iterations = int(K*std::log10(col/width));
                        if (iterations > iteration_lim) iterations = iteration_lim;
                        if (iterations < 100) iterations = 100;

                        bound.top-=height*zoom;
                        bound.bottom+=height*zoom;
                        bound.left+=width*zoom;
                        bound.right-=width*zoom;
                        generate_set();
                        output();
                }
                if (ch=='-'){
                        const long double width = bound.right-bound.left;
                        const long double height = bound.top-bound.bottom;

                        iterations = int(K*std::log10(col/width));
                        if (iterations < 100) iterations = 100;
                        if (iterations > iteration_lim) iterations = iteration_lim;

                        bound.top+=height*zoom;
                        bound.bottom-=height*zoom;
                        bound.left-=width*zoom;
                        bound.right+=width*zoom;
                        generate_set();
                        output();
                }
        }
        delete[] grid;

        refresh();
        endwin();

        return 0;}
