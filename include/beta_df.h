#ifndef BETA_DF_H
#define BETA_DF_H

// diffusion coefficient Df(dp) (ref: Kruis & Pratsinis 1993)
double Df(const double &dp);

// collision frequency function ß(dpi, dpj, regime) 
double beta(const double &dpi, const double &dpj, const int &regime);

#endif
