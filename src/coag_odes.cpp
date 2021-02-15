#include <vector>
#include "coag_odes.h"
#include "consts.h"

using namespace consts;

Coag_distr::Coag_distr( vector< vector<double> > beta_all) : 
    beta_matr(beta_all) { } // constructor

void Coag_distr::operator() ( const state_type &n , state_type &dndt, const double /* t */ ) {
    for (int k=0; k<N_ds; k++) {
        double gain = 0, loss = 0;
        for (int i=0; i<k; i++) {
            gain += beta_matr[i][k-1-i]/beta_star * n[i] * n[k-1-i]; 
            // k-1-i instead of k-i due to zero-based index
            // beta_star ß*=8kT/(3mu) just for non-dimensionalization (ref: Friedlander 2000, pg192, eqn7.18)
        }
        for (int i=0; i<N_ds; i++) {
            loss += beta_matr[i][k]/beta_star * 2 * n[k] * n[i];
        }
        dndt[k] = gain - loss;
    }
}
