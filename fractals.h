////////////////////////////////////////////////
//All fractal logic is kept here.               
//Such a cozy place.
///////////////////////////////////////////////

#ifndef _FRCTLS
#define _FRCTLS

#include <iostream>
#include <cmath>
#include <gmp.h>

enum formula{
        custom_formula,
        mandelbrot,
        julia,
        burning_ship
};

enum color{
        custom_color,
        dwell,
        DEM
};
//////////////////////////NON-PERTURBED:
template<formula f, color c>
double calculate_pixel(const double, const double);

//CUSTOM FORMULA
template<>
double calculate_pixel<custom_formula, custom_color>(const double c_r, const double c_i){
        return 0;
}

//MANDELBROT, DWELL
template<>
double calculate_pixel<mandelbrot, dwell>(const double c_r, const double c_i){
        double z_r = 0;
        double z_i = 0;
        double z2_r = 0;
        double z2_i = 0;
        double temp;
        //Periodicity
        float z_r_old = 0;
        float z_i_old = 0;
        int period = 0;
        for (int n=0; n<iters; n++){
                z2_r = z_r * z_r;
                z2_i = z_i * z_i;
                z_i = (2.0*z_i*z_r)+c_i;
                z_r = z2_r-z2_i+c_r;

                temp = std::sqrt(z2_r+z2_i);
                if (temp > bailout){
                        return n+1-log(log(temp))*inv_ln2;
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

//MANDELBROT, DEM
template<>
double calculate_pixel<mandelbrot, DEM>(const double c_r, const double c_i){
        double z_r = 0;
        double z_i = 0;
        double z2_r = 0;
        double z2_i = 0;
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
                if (temp > bailout){
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

//BURNING SHIP, DWELL
template<>
double calculate_pixel<burning_ship, dwell>(const double c_r, const double c_i){
        double z_r = 0;
        double z_i = 0;
        double z2_r = 0;
        double z2_i = 0;
        double temp;
        //Periodicity
        float z_r_old = 0;
        float z_i_old = 0;
        int period = 0;
        for (int n=0; n<iters; n++){
                z2_r = z_r * z_r;
                z2_i = z_i * z_i;
                z_i = std::abs((2.0*z_i*z_r)+c_i);
                z_r = std::abs(z2_r-z2_i+c_r);

                temp = std::sqrt(z2_r+z2_i);
                if (temp > bailout){
                        return n+1-log(log(temp))*inv_ln2;
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

//BURNING SHIP, DEM
template<>
double calculate_pixel<burning_ship, DEM>(const double c_r, const double c_i){
        double z_r = 0;
        double z_i = 0;
        double z2_r = 0;
        double z2_i = 0;
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
                z_i = std::abs((2.0*z_i*z_r)+c_i);
                z_r = std::abs(z2_r-z2_i+c_r);

                temp = std::sqrt(z2_r+z2_i);
                if (temp > bailout){
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
//////////////////////////PERTURBED:
template <formula f, color c>
double pixel_orbit(const double d_r, const double d_i, const int Z_iters);

template <>
double pixel_orbit<mandelbrot, dwell>(const double d_r, const double d_i, const int Z_iters){
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
        for (int n=0; n<Z_iters; n++){

                z2_r = z_r * z_r;
                z2_i = z_i * z_i;

                temp = z_r;
                z_r = 2*(Z_r*z_r-Z_i*z_i)+z2_r-z2_i+d_r;
                z_i = 2*(Z_r*z_i+Z_i*temp+temp*z_i)+d_i;
                temp = std::sqrt(Zz_r*Zz_r+Zz_i*Zz_i);

                if (temp > bailout){
                        return n+1-log(log(temp))*inv_ln2;
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

//MANDELBROT, DEM
template <>
double pixel_orbit<mandelbrot, DEM>(const double d_r, const double d_i, const int Z_iters){
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
        for (int n=0; n<Z_iters; n++){

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

                if (temp > bailout){
                        //extra iteration
                        dz_r_new = 2*(z_r*dz_r-z_i*dz_i)+1;
                        dz_i_new = 2*(z_r*dz_i+z_i*dz_r);
                        dz_r = dz_r_new;
                        dz_i = dz_i_new;
                        return temp*std::log(temp)/std::sqrt(dz_r*dz_r+dz_i*dz_i);
                }
                if ((z_r==z_r_old) && (z_i==z_i_old)){
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

//BURNING SHIP, DWELL
template <>
double pixel_orbit<burning_ship, dwell>(const double d_r, const double d_i, const int Z_iters){
        return 0;
}

//BURNING SHIP, DEM

template <>
double pixel_orbit<burning_ship, DEM>(const double d_r, const double d_i, const int Z_iters){
        return 0;
}

//////////////////////////REFERENCE:
int bailout2 = bailout*bailout;
template <formula f>
int reference_orbit(const mpf_t c_r, const mpf_t c_i);

//MANDELBROT
template <>
int reference_orbit<mandelbrot>(const mpf_t c_r, const mpf_t c_i){
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
                if (mpf_cmp_ui(temp, bailout2) > 0){
                        mpf_clears(z_r, z_i, z2_r, z2_i, temp, NULL);
                        return n;
                }
        }
        mpf_clears(z_r, z_i, z2_r, z2_i, temp, NULL);
        return iters;
}

//BURNING SHIP
template <>
int reference_orbit<burning_ship>(const mpf_t c_r, const mpf_t c_i){
        return 0;
}

#endif
