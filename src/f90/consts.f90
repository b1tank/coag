module consts

    ! global constants and variables 

    use preci, wp => dp
    implicit none

    real, parameter :: pi = 3.14159265358979323846
    real, parameter :: kb = 1.38064852e-23 ! Boltzmann constant
    real, parameter :: avog = 6.023e23 ! Avogadro constant
    real, parameter :: mw_air = 29.2e-3 ! molecular weight of air

    real, parameter :: mw_par = 79.866e-3 ! molecular weight of the particle material
    real, parameter :: den_par = 4230.0 ! density of the particle material

    real, parameter :: temp = 1200.0 ! temperature
    real, parameter :: visc = 1.716e-5 * ((temp/273.15)**(2.0/3.0)) ! viscosity
    real, parameter :: pres = 101325.0 ! pressure
    real, parameter :: mfp_air = visc/pres * sqrt(pi*kb*temp/(2.0*mw_air/avog)) ! mean free path of air

    real, parameter :: beta_star = 8.0*kb*temp/3.0/visc 
    ! ß* for non-dimensionalization (ref: Friedlander 2000, pg192, eqn7.18)
    integer, parameter :: regime = 0 ! regime option - default: Fuchs-Sutugin

    real, parameter :: dp1 = 0.4e-9 ! initial particle monomer diameter
    real, parameter :: n_0_dp1 = 1.0e20 ! initial monomer number concentration #/m3
    integer, parameter :: N_ds = 3000 ! ____number of concerned discrete sizes

    real, parameter :: tau = 2.0/beta_star/n_0_dp1 ! default value = 2.0837e-5 
    ! t* for non-dimensionalization, ß* instead of K !! (ref: Friedlander 2000, pg194, eqn7.24)


    ! total time is set according to table 1.1 (ref: Friedlander 2000, pg8) 
    ! n_0_dp1*t_1_10 = 1.2e16
    ! if (n_0_dp1 == 1.0e20) {
    !     t_1_10 = 1.2e16/n_0_dp1 ! = 1.2e-4
    ! }

    real, parameter :: n_star_min = 1.0/n_0_dp1 ! pyhsically miminimum total number concentration, 1 #/m3, i.e. 1e-20 after non-dimensionalization

end module consts
