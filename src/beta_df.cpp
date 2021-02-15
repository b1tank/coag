#include <iostream>
#include <cmath>
#include "consts.h"

using namespace consts;
using namespace std;


// diffusion coefficient Df(dp) (ref: Kruis & Pratsinis 1993)
double Df(const double &dp) {
    double Kn, Df;
    Kn = 2.0*lam/dp; // Knudsen number
    Df = kb*Temp/(3.0*PI*mu*dp) * ((5.0 + 4.0*Kn + 6.0*Kn*Kn + 18.0*Kn*Kn*Kn)/(5.0 - Kn + (8.0 + PI)*Kn*Kn));
    return Df;
}


// collision frequency function ß(dpi, dpj, regime) 
double beta(const double &dpi, const double &dpj, const int &regime) {
    double vi, vj, c1, c2, l1, l2, g1, g2, coeff1, coeff2, beta;
        
    vi = PI/6.0*pow(dpi, 3.0);
    vj = PI/6.0*pow(dpj, 3.0);
        
    switch(regime) {
        // Use Fuchs Sutugin interpolation expression (ref: Kruis & Pratsinis 1993)
        case(0): 
            c1 = sqrt(8.0*kb*Temp/(PI*(mw/NA)));
            c2 = c1;
            l1 = 8.0*Df(dpi)/(PI*c1);
            l2 = 8.0*Df(dpj)/(PI*c2);
            g1 = 1.0/(3.0*dpi*l1) * (pow((dpi+l1),3) - pow((dpi*dpi+l1*l1),(3.0/2))) - dpi;
            g2 = 1.0/(3.0*dpj*l2) * (pow((dpj+l2),3) - pow((dpj*dpj+l2*l2),(3.0/2))) - dpj;
            coeff1 = 2.0*PI*(dpi + dpj) * (Df(dpi) + Df(dpj));
            coeff2 = (dpi + dpj)/(dpi + dpj + 2.0*sqrt(g1*g1 + g2*g2)) + (8.0*(Df(dpi)+Df(dpj)))/((dpi+dpj)*sqrt(c1*c1+c2*c2));
            if ((dpi < dp1) || (dpj < dp1)) {
                beta = 0.0;
            } else {
                beta = coeff1 / coeff2;
            }
            break; // required and important ! ! !

        // Use Free Molecule regime (ref: Friedlander, 2000)
        case(1):
            beta = pow((3.0/(4.0*PI)),(1.0/6.0)) * pow((6.0*kb*Temp/rho_p),(1.0/2)) * pow((1.0/vi + 1.0/vj),(1.0/2)) * pow((pow(vi,(1.0/3.0)) + pow(vj,(1.0/3.0))),(2.0));
            if ((dpi < dp1) || (dpj < dp1)) {
                beta = 0.0;
            }
            break;
            
        case(2):
            beta = 1.0;
            break;

        default:
            cout << "Regime not specified" << endl;
            break;
    }
    return beta;
}
