#ifndef CONSTS_H
#define CONSTS_H

// global constants and variables (only forward declarations !!!)
namespace consts {
    extern const double PI;
    extern const double kb; // Boltzmann constant
    extern const double NA; // Avogadro constant
    extern const double mw; // molecular weight of the particle material
    extern const double mwg; // molecular weight of air
    extern const double rho_p; // density of the particle material
    extern const double Temp; // temperature
    extern const double mu; // viscosity
    extern const double P; // pressure
    extern const double lam; // mean free path of air

    extern const double beta_star; 
    // ß* for non-dimensionalization (ref: Friedlander 2000, pg192, eqn7.18)
    extern const int regime; // regime option

    extern const double dp1; // initial particle monomer diameter
    extern const double n_0_dp1; // initial monomer number concentration
    extern const int N_ds; // nunmber of concerned discrete sizes

    extern const double tau; 
    // t* for non-dimensionalization, ß* instead of K !! (ref: Friedlander 2000, pg194, eqn7.24)

    extern const double n_star_min;
}

#endif
