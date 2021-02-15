/* consts.cpp */
#include <cmath>
// global constants and variables 
// only instantiated once to save rebuild time !
// however, these constants are no longer compile-time available !!

namespace consts {
    extern const double PI = 3.14159265358979323846;
    extern const double kb = 1.38064852e-23; // Boltzmann constant
    extern const double NA = 6.023e23; // Avogadro constant
    extern const double mwg = 29.2e-3; // molecular weight of air

    extern const double mw = 79.866e-3; // molecular weight of the particle material
    extern const double rho_p = 4230.0; // density of the particle material

    extern const double Temp = 1200.0; // temperature
    extern const double mu = 1.716e-5 * pow((Temp/273.15),(2.0/3.0)); // viscosity
    extern const double P = 101325.0; // pressure
    extern const double lam = mu/P * sqrt(PI*kb*Temp/(2.0*mwg/NA)); // mean free path of air

    extern const double beta_star = 8.0*kb*Temp/3.0/mu; 
    // ß* for non-dimensionalization (ref: Friedlander 2000, pg192, eqn7.18)
    extern const int regime = 0; // regime option - default: Fuchs-Sutugin

    extern const double dp1 = 0.4e-9; // initial particle monomer diameter
    extern const double n_0_dp1 = 1.0e20; // initial monomer number concentration #/m3
    extern const int N_ds = 20000; // ____nunmber of concerned discrete sizes

    extern const double tau = 2.0/beta_star/n_0_dp1; // default value = 2.0837e-5 
    // t* for non-dimensionalization, ß* instead of K !! (ref: Friedlander 2000, pg194, eqn7.24)


    /* total time is set according to table 1.1 (ref: Friedlander 2000, pg8) 
    n_0_dp1*t_1_10 = 1.2e16
    if (n_0_dp1 == 1.0e20) {
        t_1_10 = 1.2e16/n_0_dp1; // = 1.2e-4
    }
    */

    extern const double n_star_min = 1.0/n_0_dp1; // pyhsically miminimum total number concentration, 1 #/m3, i.e. 1e-20 after non-dimensionalization
}
