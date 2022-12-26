//
//  MGLM.hpp
//
//
//  Created by Kuo-Jung Lee on 2022/3/17.
//
//#include <math.h>

#include <stdio.h>
//#include <algorithm>
//#include <assert.h>
#include <cmath>
#include <ctime>    // For time()
//#include <cstdlib>  // For srand() and rand()
//#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <iomanip>
//#include <list>
//#include <limits>
#include <vector>
#include <string>
//#include <sstream>
#include <algorithm>


#include<iostream>
//#include<chrono> //for sleeping
#include<thread> // --do--
#include<cstdlib>//for random increments


//#include <RcppArmadillo.h>
#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

using namespace std;

using namespace Rcpp;
using namespace arma;


#ifndef MGLM_hpp
#define MGLM_hpp

#include <stdio.h>

#include "tmvrnormGibbs_KJLEE.h"
class MultiLinearModel{
//private:
    //mat Ri0, Ri1, Ri_all, Ri0_inv, Ri1_inv, Ri_all_inv;
protected:
    int Num_of_iterations;
    int Num_of_obs, Num_of_FixedEffs, Num_of_RandomEffs, Num_of_attributes;
    int Num_of_deltas, Num_of_alphas;
    cube X, WW, VV, Z;
    mat Y, Y_Binary ;
    
    field<cube> W, V; //(Nxa)xTxT
    
    vec TimePointsAvailable;//group_indices,

    mat beta_samples, alpha_samples, delta_samples;
    cube b_samples, Sigma_samples;
    
    double sigma2_beta, Vb, sigma2_alpha, sigma2_delta;
    
    double tuning_alpha, tuning_delta;
    mat Lambda; 
    
    List Data, InitialValues, HyperPara, UpdatePara, TuningPara;
    
    vec beta_mean, alpha_mean, delta_mean;
    mat Sigma_mean, b_mean;
    
    double PCP, PRE, ePCP, logL, AIC, cAIC, BIC, DIC, MPL; 
    
    bool updateystar, updatebeta, updateb, updateSigma, updatealpha, updatedelta;
    
    double acc_rate_alpha, acc_rate_delta; 
    
public:
    MultiLinearModel(int iNum_of_iterations, List list_Data, List list_InitialValues, List list_HyperPara, List list_UpdatePara, List list_TuningPara);
    
    //mat Ri(int tp, cube x, vec delta);
    mat Ri_Within(int i, int tp,  vec delta);
    mat Ri_Cross(int i, int tp,  vec delta);
    
    void Update_ystar(int iter);
    void Update_b(int iter);
    void Update_beta(int iter);
    void Update_Sigma(int iter);
    void Update_alpha(int iter);
    void Update_delta(int iter);
    void Update_delta_CompWise(int iter); 
    void ParameterEstimation(); 
    SEXP MCMC_Procedure();
};



#endif /* MGLM_hpp */
