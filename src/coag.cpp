#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <time.h>
#include <cmath>
#include <boost/numeric/odeint.hpp>
#include "consts.h"
#include "beta_df.h"
#include "coag_odes.h"

using namespace std;
using namespace consts;
using namespace boost::numeric::odeint;


/********************************* Main Program -- Start *********************************/
int main() {

    /******************** Clock -- Start ********************/
    clock_t t_cpu; 
    t_cpu = clock(); // cpu time observer

    /******************** Reading input  parameters ********************/
    ifstream infile;
    infile.open("../input.txt");

    //double tstep, tf;
    double tstep;
    //bool observer;
    //infile >> observer;
    infile >> tstep;
    //infile >> tf;

    infile.close();

    /******************** Initial value conditions ********************/
    // number concentration initial values at t = 0
    state_type n_k;
    n_k.push_back(1.0); // n_k[0]=1.0 after being non-dimensionalized by n_0_dp1
    for (int i=1; i<N_ds; i++) {
        n_k.push_back(0.0);
    }

    /******************** Parameters initialization ********************/
    // all dp can be initialized once N_ds is given (dimer, trimer, tetramer ...)
    vector<double> dp;
    vector<int> dp_out;
    dp.push_back(dp1); // dp[0] = dp1 -> watch out the index trick !
    double v1 = PI/6.0*dp1*dp1*dp1;
    for (int i=1; i<N_ds; i++) {
        dp.push_back(pow(v1*(i+1)*6.0/PI, 1.0/3.0));
    }
    for (int i=1000; i<=N_ds; i += 1000) {
        dp_out.push_back(i);
    }

    // all beta can be initialized (permanently) in a matrix once all dp are known
    vector< vector<double> > beta_matr;
    for (int i=0; i<N_ds; i++) {
        vector<double> beta_tmp;
        for (int j=0; j<N_ds; j++) {
            beta_tmp.push_back(beta(dp[i],dp[j],regime));
        }
        beta_matr.push_back(beta_tmp);
    }

    /* print out all ß(dpi,dpj) for debugging
    for (int i=0; i<N_ds; i++) {
        for (int j=0; j<N_ds; j++) {
            cout << beta_matr[i][j] << ' '; 
        }
        cout << endl;
    }
    */


    /******************** Defining and solving ODEs ********************/
    // define ODEs using all ß(dpi,dpj)
    Coag_distr coag(beta_matr);

    // observer vectors for n_k and time
    vector< state_type > n_k_vec;
    vector< double > times;

    ///******************** integrate function with adaptive step size ********************/
    //size_t steps;
    //// solve ODEs use integrate function (adaptive step size) in boost library
    //if (observer) {
    //    steps = integrate( coag, n_k, 0.0, tf/tau, tstep/tau, push_back_state_and_time( n_k_vec , times ) );
    //}else {
    //    steps = integrate( coag, n_k, 0.0, tf/tau, tstep/tau );
    //}
    ///******************** single-stepper with constant step size ********************/
    // runge_kutta_cash_karp54< state_type > rkstepper; // constant step size
    size_t steps = 0;
    double t = 0.0;
    double st = tstep/tau;
    //******************** single-stepper with controlled (adaptive) step size ********************/
    typedef runge_kutta_cash_karp54< state_type > error_stepper_type; 
    typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type; // customize controlled stepper using error stepper
    double abs_err = 1.0e-6 , rel_err = 1.0e-6 , a_x = 1.0 , a_dxdt = 1.0;
    controlled_stepper_type controlled_stepper( default_error_checker< double , range_algebra , default_operations >( abs_err , rel_err , a_x , a_dxdt ) );

    // initialize observer
    push_back_state_and_time observer( n_k_vec , times );
    observer( n_k , t ); 

    int j = 0; // output iterator
    int count; // counting number of trailing "0"s in a certain n_k vector
    int flag_ds; // exact number of discrete sizes where particles start emerging
    while (n_k[N_ds-1] <= n_star_min) {
        // rkstepper.do_step( coag, n_k, t, st ); // constant step size
        if (controlled_stepper.try_step( coag, n_k, t, st )) { // "1" means step size rejected while "0" means accepted ! ! !
            //cout << "fail: " << t << ' ' << st << endl;
            continue; // if step size is rejected, try smaller step size
        }
        //cout << "succ: " << t << ' ' << st << endl;
        observer( n_k , t ); 
        steps += 1;
        if (n_k[dp_out[j]-1] >= n_star_min ) {
            count = 0;
            for (int k=N_ds; k>=0; k--) {
                if (n_k_vec[steps-1][k] >= n_star_min) {
                    break;
                }
                count += 1;
            }
            flag_ds = N_ds - count - 1;
            string filename = "../output/output" + to_string(flag_ds) + ".txt";
            ofstream outfile (filename);
            if (outfile.is_open()) {
                outfile << "0" << '\t';
                for (int i=0; i<dp_out[j]; i++) {
                    outfile << dp[i] * 1.0e9 << '\t';
                }
                outfile << endl;
                outfile << times[steps-1]*tau << '\t';
                for (int k=0; k<dp_out[j]; k++) {
                    if (n_k_vec[steps-1][k] < n_star_min) {
                        outfile << 0.0 << '\t'; // non-dimensionalized
                    }else {
                        outfile << n_k_vec[steps-1][k] << '\t'; // non-dimensionalized
                    }
                }
                outfile << endl;
                outfile.close();
                cout << "output" << flag_ds << ".txt is ready..."; 
                printf("CPU time: %f s; ", ((float)(clock()-t_cpu))/CLOCKS_PER_SEC);
                printf("Coag time: %10.8f s\n", t*tau);
            }
            j++;
        }
    }

      /******************** Output to file ********************/
    ofstream outfile ("../output/output.txt");
    if (outfile.is_open()) {
        outfile << "0" << '\t';
        for (int i=0; i<N_ds; i++) {
            outfile << dp[i] * 1.0e9 << '\t';
        }
        outfile << endl;

        //if (observer) {
        //    cout << "Observer on!" << endl;
            for (size_t i=0; i<=steps; i++) {
                outfile << times[i] * tau << '\t'; // dimensionalizing 
                for (int j=0; j<N_ds; j++) {
                    if (n_k_vec[i][j] < n_star_min) {
                        outfile << 0.0 << '\t'; // non-dimensionalized
                    }else {
                        outfile << n_k_vec[i][j] << '\t'; // non-dimensionalized
                    }
                    //outfile << n_k_vec[i][j] * n_0_dp1 << '\t'; // dimensionalizing
                }
                outfile << endl;
            }
        //}else {
        //    cout << "Observer off!" << endl;
        //    for (int i=0; i<N_ds; i++) {
        //            outfile << n_k[i] << '\t'; // non-dimensionalized
        //    }
        //    outfile << endl;
        //}

        // check the total number concentration at the final time
        double sn = 0; 
        for (int j=0; j<N_ds; j++) {
            sn += n_k_vec[steps][j];
        }
        cout << "n_tot/n_0_tot: " << sn/1.0 << endl; 
        cout << "Coag time: " << t << " (non dimensional)" << endl; // = 12, reach self-preserving (Hidy, 1965)
        printf("Coag time: %10.8f s\n", t*tau);
        // check the total number concentration at the final time

        outfile.close();
    }

    /******************** Clock -- End ********************/
    t_cpu = clock() - t_cpu;
    printf("CPU time: %f s\n", ((float)t_cpu)/CLOCKS_PER_SEC);

    return 0;
}
/********************************* Main Program -- End *********************************/
