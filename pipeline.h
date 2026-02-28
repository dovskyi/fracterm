///////////////////////////////////////
//Too big to put in main, cool mariani
//silver shenanigans
///////////////////////////////////////
#ifndef PIPE
#define PIPE

#include <iostream>
#include <cmath>

template <bool perturb, formula f, color c>
char process_pixel(const double& dy, const double& dx, const double& c_real, const double& c_imag, const perturb_data* const p_d){
        double modulus;
        int cnt = 0;
        double mu = 0;
        //antialiasing yaay
        for (int i=0; i<2; i++){
                double offset_real = c_real+offset[i]*dy;
                for (int j=0; j<1; j++){
                        double offset_imag = c_imag+offset[j]*dx;
                        if constexpr(perturb){
                                //technically c here would be delta
                                modulus = pixel_orbit<f, c>(offset_real, offset_imag, p_d->Z_iters);
                        }
                        else {
                                modulus = calculate_pixel<f, c>(offset_real, offset_imag);
                        }
                        if (modulus > 0){
                                cnt++;
                                mu+=modulus;
                        }
                }
        }
        if (cnt > 1){
                if constexpr(c == DEM){
                        double norm = (mu*glow)/(cnt*dx);
                        if (norm > norm_lim) norm = norm_lim;
                        double index = std::log(1.0+norm)*inv_ln;
                        //this +0.1 works on hopes and dreams:
                        //i hope and dream it does not segfault
                        return density[int((1.0-index)*(d_size+0.1))];
                        //it probably mathematically can't? idk
                }
                if constexpr(c == dwell){
                        double norm = (mu/cnt)*0.07;
                        return density[int(norm)%(d_size+1)];
                }
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


template <bool perturb, formula f, color c>
void initial_perimeter(const double& dy, const double& dx, const mariani_data& m_d, const perturb_data* const p_d){
        if constexpr(perturb){
                double d_real;
                double d_imag;
                double temp;
                for (int x=0; x<row; x++){
                        temp = x-p_d->iZ_i;
                        d_imag = temp*dx;

                        temp = -p_d->iZ_r;
                        d_real = temp*dy;
                        grid[x*col] = process_pixel<perturb, f, c>(dy, dx, d_real, d_imag, p_d);
                        temp = p_d->iZ_r;
                        d_real = temp*dy;
                        grid[x*col+col-1] = process_pixel<perturb, f, c>(dy, dx, d_real, d_imag, p_d);
                }
                for (int y=0; y<col; y++){
                        temp = y-p_d->iZ_r;
                        d_real= temp*dy;

                        temp = -p_d->iZ_i;
                        d_imag = temp*dx;
                        grid[y] = process_pixel<perturb, f, c>(dy, dx, d_real, d_imag, p_d);
                        temp = p_d->iZ_i;
                        d_imag = temp*dx;
                        grid[(row-1)*col+y] = process_pixel<perturb, f, c>(dy, dx, d_real, d_imag, p_d);
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
                        grid[x*col] = process_pixel<perturb, f, c>(dy, dx, c_real, c_imag, p_d);
                        c_real = left_d+((col-1)*dy);
                        grid[x*col+col-1] = process_pixel<perturb, f, c>(dy, dx, c_real, c_imag, p_d);
                }

                for (int y=0; y<col; y++){
                        c_real = left_d+(y*dy);

                        c_imag = bottom_d+dx;
                        grid[y] = process_pixel<perturb, f, c>(dy, dx, c_real, c_imag, p_d);
                        c_imag = bottom_d+((row-1)*dx);
                        grid[(row-1)*col+y] = process_pixel<perturb, f, c>(dy, dx, c_real, c_imag, p_d);
                }
        }
}

template <bool perturb, formula f, color c>
void mariani_silver(const double& dy, const double& dx, const mariani_data m_d, const perturb_data* const p_d){
        const int round_h = m_d.height/2;
        const int round_w = m_d.width/2;
        const int midp_h = m_d.height-round_h;
        const int midp_w = m_d.width-round_w;
        if constexpr(f==burning_ship){
        //small amendement for burning ship, lower part fills otherwise
        if (check_perimeter(m_d) == 1 && dx < 0.0005) {
                for (int y=0; y<m_d.width-0; y++){
                        for (int x=0; x<m_d.height-0; x++){
                                grid[(m_d.h_offset+x)*col+(m_d.w_offset+y)] = grid[m_d.h_offset*col+m_d.w_offset];
                        }
                }
                return;
        }
        }
        else{
        if (check_perimeter(m_d) == 1 && dx < 0.01) {
                for (int y=0; y<m_d.width-0; y++){
                        for (int x=0; x<m_d.height-0; x++){
                                grid[(m_d.h_offset+x)*col+(m_d.w_offset+y)] = grid[m_d.h_offset*col+m_d.w_offset];
                        }
                }
                return;
        }
        }

        if constexpr(perturb){
                //d is for delta btw (btw)
                double d_real;
                double d_imag;
                double temp;
                if (m_d.height <= 4 || m_d.width <= 4){
                        for (int y=0; y<m_d.width; y++){
                                temp = (y+m_d.w_offset)-p_d->iZ_r;
                                d_real = temp*dy;
                                for (int x=0; x<m_d.height;x++){
                                        temp = (x+m_d.h_offset)-p_d->iZ_i;
                                        d_imag = temp*dx;
                                        grid[(m_d.h_offset+x)*col+(m_d.w_offset+y)] = process_pixel<perturb, f, c>(dy, dx, d_real, d_imag, p_d);
                                }
                        }
                        return;
                }

                temp = (m_d.w_offset+round_w)-p_d->iZ_r;
                d_real = temp*dy;
                for (int x=0; x<m_d.height; x++){
                        temp = (x+m_d.h_offset)-p_d->iZ_i;
                        d_imag = temp*dx;
                        grid[(m_d.h_offset+x)*col+(m_d.w_offset+round_w)] = process_pixel<perturb, f, c>(dy, dx, d_real, d_imag, p_d);
                }

                temp = (m_d.h_offset+round_h)-p_d->iZ_i;
                d_imag = temp*dx;
                for (int y=0; y<m_d.width; y++){
                        temp = (y+m_d.w_offset)-p_d->iZ_r;
                        d_real = temp*dy;
                        grid[(m_d.h_offset+round_h)*col+(m_d.w_offset+y)] = process_pixel<perturb, f, c>(dy, dx, d_real, d_imag, p_d);
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
                                        grid[(m_d.h_offset+x)*col+(m_d.w_offset+y)] = process_pixel<perturb, f, c>(dy, dx, c_real, c_imag, p_d);
                                }
                        }
                        return;
                }
                c_real = left_d+((m_d.w_offset+round_w)*dy);
                for (int x=0; x<m_d.height; x++){
                        c_imag = bottom_d+((m_d.h_offset+x)*dx);
                        grid[(m_d.h_offset+x)*col+(m_d.w_offset+round_w)] = process_pixel<perturb, f, c>(dy, dx, c_real, c_imag, p_d);
                }
                c_imag = bottom_d+((m_d.h_offset+round_h)*dx);
                for (int y=0; y<m_d.width; y++){
                        c_real = left_d+((m_d.w_offset+y)*dy);
                        grid[(m_d.h_offset+round_h)*col+(m_d.w_offset+y)] = process_pixel<perturb, f, c>(dy, dx, c_real, c_imag, p_d);
                }
        }

        //call 4 new rectangles
        mariani_silver<perturb, f, c>(dy, dx, {m_d.h_offset        , m_d.w_offset        , round_h, round_w}, p_d);
        mariani_silver<perturb, f, c>(dy, dx, {m_d.h_offset        , m_d.w_offset+round_w, round_h, midp_w }, p_d);
        mariani_silver<perturb, f, c>(dy, dx, {m_d.h_offset+round_h, m_d.w_offset        , midp_h , round_w}, p_d);
        mariani_silver<perturb, f, c>(dy, dx, {m_d.h_offset+round_h, m_d.w_offset+round_w, midp_h , midp_w }, p_d);
        return;
}
#endif
