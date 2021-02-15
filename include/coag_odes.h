#ifndef COAG_ODES_H
#define COAG_ODES_H

#include <vector>
using namespace std;

typedef vector<double> state_type; // storing calculated values for the vector of variables in the ode solver

// Intermediate timesteps observer
struct push_back_state_and_time {
    vector< state_type >& m_states;
    vector< double >& m_times;

    push_back_state_and_time( vector< state_type >& states, vector< double >& times ) : m_states( states ) , m_times( times ) { }

    void operator()( const state_type &x, double t ) {
        m_states.push_back( x );
        m_times.push_back( t );
    }
};


// rhs of dynamic equation of discrete coagulation-only distribution (non-dimensionalized)
class Coag_distr {

    vector< vector<double> > beta_matr; // N_ds * N_ds matrix of ß values calculated beforehand

public:
    Coag_distr( vector< vector<double> > beta_all); // constructor
    void operator() ( const state_type &n , state_type &dndt, const double /* t */ ); 
};

#endif
