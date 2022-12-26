//
//  MGLM.cpp
//
//
//  Created by kuojung on 2022/3/17.
//

// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "MGLM.h"
//RNGScope scope;
//#include <RcppArmadilloExtensions/sample.h>

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

MultiLinearModel::MultiLinearModel(int iNum_of_iterations, List list_Data, List list_InitialValues, List list_HyperPara, List list_UpdatePara, List list_TuningPara)
{
    Num_of_iterations = iNum_of_iterations;
    Data = list_Data;
    InitialValues = list_InitialValues;
    HyperPara = list_HyperPara;
    UpdatePara = list_UpdatePara;
    TuningPara = list_TuningPara;
    
    updateystar = as<bool>(UpdatePara["UpdateYstar"]);
    updatebeta = as<bool>(UpdatePara["UpdateBeta"]);
    updateb = as<bool>(UpdatePara["Updateb"]);
    updateSigma = as<bool>(UpdatePara["UpdateSigma"]);
    updatedelta = as<bool>(UpdatePara["UpdateDelta"]);
    updatealpha = as<bool>(UpdatePara["UpdateAlpha"]);

    

    
    acc_rate_alpha = 0;
    acc_rate_delta = 0;
    
    Rcout<< "Read Data" << endl;
    //Y = as<mat>(Data["Y"]); //N*(TxK)
    
    Y_Binary = as<mat>(Data["Y"]); //N*(TxK)
    
    Y = Y_Binary;
    
    //Rcout << "size(Y_Binary) ="  <<size(Y_Binary) << endl;
    //Rcout << "Y_binary =" << endl << Y_Binary.submat(0, 0, 9, 9) << endl;
    
    //Rcout << size(Y) << endl;
    X = as<cube>(Data["X"]); // N x (TxP) x (KxP)
    Z = as<cube>(Data["Z"]);
    
    
    
    //Rcout << "size(X) = "<<size(X) << endl;
    
    //mat X_tmp = X.subcube(0, 0, 0, 19, 0, 39);
    
    //Rcout << X_tmp << endl;

    
    
    //Rcout << "size(Z) = " << size(Z) << endl;
    
    //mat Z_tmp = Z.subcube(0, 0, 0, 0, 19, 4);
    
    //Rcout << Z_tmp << endl;

    TimePointsAvailable = as<vec>(Data["TimePointsAvailable"]);
    //X_all = as<cube>(Data["X.all"]);
    //group_indices = as<vec>(Data["group.indices"]);
    
    Num_of_obs = Y.n_cols;
    Num_of_attributes = X.n_cols/TimePointsAvailable.max();
    Num_of_FixedEffs = X.n_slices/Num_of_attributes;
    Num_of_RandomEffs = Z.n_slices/Num_of_attributes;
    
    Num_of_alphas = 0;
    Num_of_deltas = 0;
    
    if(updatealpha){
        //Rcout << "Check 1" << endl;
        W.reset();
        WW.reset();
        WW = as<cube>(Data["W"]);
        Num_of_alphas = WW.n_rows/Num_of_obs;
        W.set_size(Num_of_obs);
        for(int i = 0; i<Num_of_obs; i++)
            W(i) = WW.rows( (i*Num_of_alphas), (Num_of_alphas*(i+1)-1));
        alpha_samples.set_size(Num_of_alphas, Num_of_iterations);
        alpha_samples.zeros();
        alpha_mean.zeros(Num_of_alphas);
        alpha_samples.col(0) = as<vec>(InitialValues["alpha.ini"]);
        sigma2_alpha = as<double>(HyperPara["sigma2.alpha"]);
        
        //Rcout << "size(WW) = " << size(WW) << endl;
        //Rcout<< "Read Tuning parameters." << endl;
        tuning_alpha = as<double>(TuningPara["tuning.alpha"]);
        Rcout << "W(i) = " << endl << size(W(0)) << endl;

    }
    
    
    
    if(updatedelta){
        V.reset();
        VV.reset();
        VV = as<cube>(Data["V"]);
        Num_of_deltas = VV.n_rows/Num_of_obs;
        V.set_size(Num_of_obs);
        for(int i = 0; i<Num_of_obs; i++){
            V(i) = VV.rows( (i*Num_of_deltas), (Num_of_deltas*(i+1)-1));
        }

        delta_samples.set_size(Num_of_deltas, Num_of_iterations);
        delta_samples.zeros();
        delta_mean.zeros(Num_of_deltas);
        delta_samples.col(0) = as<vec>(InitialValues["delta.ini"]);
        sigma2_delta = as<double>(HyperPara["sigma2.delta"]);
        //Rcout << "size(VV) = " << size(VV) << endl;
        
        //Rcout<< "Read Tuning parameters." << endl;
        
        tuning_delta = as<double>(TuningPara["tuning.delta"]);
        Rcout << "V(i) = " << endl << size(V(0)) << endl;
    }
    
    
    


//resize the posterior sample size
    
    //y_star.set_size( size(Y) );
    
    //Rcout << "size of y.star = " << size(y_star) << endl;
    
    
    beta_samples.set_size((Num_of_FixedEffs*Num_of_attributes), Num_of_iterations);
    beta_samples.zeros();
    beta_mean.zeros((Num_of_FixedEffs*Num_of_attributes));
    

    Sigma_samples.set_size(Num_of_attributes*Num_of_RandomEffs, Num_of_attributes*Num_of_RandomEffs, Num_of_iterations);
    Sigma_samples.slice(0).eye();
    
    Sigma_mean.eye(Num_of_attributes*Num_of_RandomEffs, Num_of_attributes*Num_of_RandomEffs);
    
    b_samples.set_size(Num_of_obs, Num_of_attributes*Num_of_RandomEffs, Num_of_iterations);
    b_samples.zeros();
    b_mean.zeros(Num_of_obs, Num_of_attributes*Num_of_RandomEffs); 
    
//===========================================================================//
    
    Rcout<< "Initial Values" << endl;
    beta_samples.col(0) = as<vec>(InitialValues["beta.ini"]);
    b_samples.slice(0) = as<mat>(InitialValues["b.ini"]);
    Sigma_samples.slice(0) = as<mat>(InitialValues["Sigma.ini"]);

    

  
    Rcout<< "Read Hyperparameters." << endl;
    // Hyperparameters
    sigma2_beta = as<double>(HyperPara["sigma2.beta"]);
    Vb = as<double>(HyperPara["Vb"]);
    Lambda = as<mat>(HyperPara["Lambda"]);




    
    
    
    Rcout<< "All set" << endl;
}


/*
mat MultiLinearModel::Ri(int tp, cube x, vec delta)
{
    //Rcout << "delta = " << delta << endl;
    mat F_tmp(tp, tp), F(tp, tp);
    F.zeros();
    F_tmp.zeros();
    vec x_vec;
    
    //Rcout << "tp = " << tp << endl;
    //Rcout << "size(x)= " << size(x) << endl;

    for(int i = 0; i<tp; i++)
        for(int j=0; j<tp; j++){
            //Rcout << "i = " << i << "\tj = " << j << endl;
            x_vec = x.slice(i).col(j);
            //Rcout << "x_vec = " << x_vec << endl;
            //Rcout << "delta = " << delta << endl;
            F_tmp(i, j) = accu(delta.t()*x_vec);
        }
            
    F_tmp = datum::pi*exp(F_tmp)/(1.+exp(F_tmp));
    F(0, 0) = 1;
    
    
    //Rcout << "F_tmp = \n" << F_tmp << endl;
    
    for(int t=1; t<tp; t++)
        F(t, 0) = cos(F_tmp(t, 0));
    for(int j = 1; j<tp-1; j++)
        for(int t = j+1; t<tp; t++)
            F(t, j) = cos(F_tmp(t, j))*prod(sin(F_tmp(t, span(0, j-1) )));
    for(int t=1; t<tp; t++)
        F(t, t) = prod(sin(F_tmp(t, span(0, t-1) )));
    mat Ri = F * F.t();
    
    //Rcout << "Done Ri" << endl;
    return (Ri);
}
*/

mat MultiLinearModel::Ri_Within(int i, int tp,  vec alpha)
{
    //Rcout << "alpha = \n" << alpha << endl;
    mat F_tmp(tp, tp), F(tp, tp);
    F.zeros();
    F_tmp.zeros();
    
    mat  W_tmp;

    for(int delta_ind_U = 0; delta_ind_U<alpha.n_elem; delta_ind_U++){
        W_tmp = W(i).row(delta_ind_U);
        //Rcout << "delta_ind_U = " << delta_ind_U << endl;
        //Rcout << "as_scalar(alpha(delta_ind_U))*W_tmp = \n" << as_scalar(alpha(delta_ind_U))*W_tmp << endl;
        F_tmp += as_scalar(alpha(delta_ind_U))*W_tmp(0, 0, size(tp, tp) ); //W(delta_ind_U).row(i); //slice(i)(0, 0, size(tp, tp) );
        //Rcout << "Done within" << endl;
    }
    
    F_tmp = datum::pi*exp(F_tmp)/(1.+exp(F_tmp));
    F(0, 0) = 1;
    
    for(int t=1; t<tp; t++)
        F(t, 0) = cos(F_tmp(t, 0));
    for(int j = 1; j<tp-1; j++)
        for(int t = j+1; t<tp; t++)
            F(t, j) = cos(F_tmp(t, j))*prod(sin(F_tmp(t, span(0, j-1) )));
    for(int t=1; t<tp; t++)
        F(t, t) = prod(sin(F_tmp(t, span(0, t-1) )));
    mat Ri = F * F.t();
    return (Ri);
}

mat MultiLinearModel::Ri_Cross(int i, int tp,  vec delta)
{
    //Rcout << "delta = \n" << delta << endl;
    mat F_tmp(tp, tp), F(tp, tp);
    F.zeros();
    F_tmp.zeros();
    
    mat V_tmp;

    for(int delta_ind_U = 0; delta_ind_U<delta.n_elem; delta_ind_U++){
        V_tmp = V(i).row(delta_ind_U);
        //Rcout << "delta_ind_U = " << delta_ind_U << endl;
        //Rcout << "as_scalar(delta(delta_ind_U))*V_tmp = \n" << as_scalar(delta(delta_ind_U))*V_tmp << endl;
        F_tmp += as_scalar(delta(delta_ind_U))*V_tmp(0, 0, size(tp, tp) ); //V(delta_ind_U).slice(i)(0, 0, size(tp, tp) );
        //Rcout << "Done Cross " << endl;
    }
    F_tmp = datum::pi*exp(F_tmp)/(1.+exp(F_tmp));
    F(0, 0) = 1;
    
    for(int t=1; t<tp; t++)
        F(t, 0) = cos(F_tmp(t, 0));
    for(int j = 1; j<tp-1; j++)
        for(int t = j+1; t<tp; t++)
            F(t, j) = cos(F_tmp(t, j))*prod(sin(F_tmp(t, span(0, j-1) )));
    for(int t=1; t<tp; t++)
        F(t, t) = prod(sin(F_tmp(t, span(0, t-1) )));
    mat Ri = F * F.t();
    return (Ri);
}


void MultiLinearModel::Update_ystar(int iter)
{
    //Rcout << "Update ystar" << endl;
    mat Ri1, Ri0, Ri_all;
    mat X_tmp, Z_tmp;
    vec y_mu;
    int tp, tp_all;
    vec lower, upper;
    vec y_tmp;
    for(int i=0; i<Num_of_obs; i++){
        //Rcout << "i = " << i << endl;
        tp = TimePointsAvailable(i);
        tp_all = tp*Num_of_attributes;

        if(updatealpha)
            Ri1 = Ri_Within(i, TimePointsAvailable(i), alpha_samples.col(iter));
        else
            Ri1.eye(TimePointsAvailable(i), TimePointsAvailable(i));
        
        if(updatedelta)
            Ri0 = Ri_Cross(i, Num_of_attributes, delta_samples.col(iter));
        else
            Ri0.eye(Num_of_attributes, Num_of_attributes);

        
        //Ri1 = Ri(tp, W(i), alpha_samples.col(iter));
        //Ri0 = Ri(Num_of_attributes, V(i), delta_samples.col(iter));
        Ri_all = kron(Ri1, Ri0);
        
        //Rcout << "Check 0" << endl;
        X_tmp = X(span(i), span(0, tp_all-1), span::all); //X.row(i);
        Z_tmp = Z(span(i), span(0, tp_all-1), span::all); //Z.row(i);
                
        //Rcout << "Check 1" << endl;
        
        //lower.elem(find(Y(span(0, tp-1), i)>0)).zeros();
        //lower.elem(find(Y(span(0, tp-1), i)==0)).ones();
        //lower.elem(find(Y(span(0, tp-1), i)==0)) *= -datum::inf;
        
        //upper.elem(find(Y(span(0, tp-1), i)==0)).zeros();
        //upper.elem(find(Y(span(0, tp-1), i)>0)).ones();
        //upper.elem(find(Y(span(0, tp-1), i)>0)) *= datum::inf;
        
        //Rcout << "check 1" << endl;
        //X_tmp = X(span(0, tp-1), span(0, Num_of_covariates-1), span(i));
        //Z_tmp = (Z.slice(i).rows(0, tp-1));
        //b_tmp = vectorise(b_samples(i, span(0, tp_all-1), iter));
        //Rcout << "beta_samples.col(iter) = " << beta_samples.col(iter) << endl;
        
        //Rcout << "b_tmp = " << b_tmp << endl;
        
        //Rcout << " X_tmp* beta_samples = " << endl << X_tmp* beta_samples.col(iter) << endl;
        //Rcout << "Z_tmp*b_tmp = " << endl << Z_tmp*b_tmp << endl;
        y_mu = X_tmp* beta_samples.col(iter)+Z_tmp*b_samples.slice(iter).row(i).t();
            
            //.slice(iter).row(i).t();;
        
        //Rcout << "check 2" << endl;

        lower.set_size(tp_all);
        upper.set_size(tp_all);
        
        

        lower.elem(find(Y_Binary(span(0, tp_all-1), i)>0)).zeros();
        lower.elem(find(Y_Binary(span(0, tp_all-1), i)==0)).ones();
        lower.elem(find(Y_Binary(span(0, tp_all-1), i)==0)) *= -datum::inf;
        
        upper.elem(find(Y_Binary(span(0, tp_all-1), i)==0)).zeros();
        upper.elem(find(Y_Binary(span(0, tp_all-1), i)>0)).ones();
        upper.elem(find(Y_Binary(span(0, tp_all-1), i)>0)) *= datum::inf;

        //Rcout << "check 3" << endl;
        
        //Rcout << "y_mu =" << size(y_mu) << endl;
        //Rcout << "Ri_all =" << size(Ri_all) << endl;
        //Rcout << "Sample = " << rtmvnorm_gibbs_KJLEE(1, y_mu, Ri_all, lower, upper, 100, zeros<vec>(tp), 5).t() << endl;
        //Rcout << size(Y) << endl;
        
        Y(span(0, tp_all-1), i) = rtmvnorm_gibbs_KJLEE(1, y_mu, Ri_all, lower, upper, 100, zeros<vec>(tp*Num_of_attributes), 5).t();
        
        //Rcout << "check 4" << endl;
        if(Y.col(i).has_nan()){
            y_tmp = Y.col(i);
            y_tmp(find_nonfinite(y_tmp)).zeros();
            Y.col(i) = y_tmp;
        }
        //Rcout << "check 4" << endl;
 
        Y.col(i) = clamp(Y.col(i), -10, 10);
    }

    //Rcout << "Done Update ystar" << endl;
}

void MultiLinearModel::Update_beta(int iter)
{
    //Rcout<< "Generate beta" <<endl;

    mat Ri1_tmp, Ri0_tmp, Ri_all_inv;
    
    int beta_size = Num_of_attributes*Num_of_FixedEffs;
    mat beta_cov= zeros<mat>(beta_size,beta_size);
    vec beta_mu = zeros<vec>(beta_size);
    mat X_tmp, Z_tmp;
    vec Res;
    vec tmp;
    
    int tp, tp_all;
    
    mat Ind_mat_beta = eye(Num_of_attributes*Num_of_FixedEffs, Num_of_attributes*Num_of_FixedEffs);
    //beta_cov =X.col(0);
    //Rcout << b_samples.slice(0).row(0) << endl;
    for(int i=0; i<Num_of_obs; i++){
        
        tp = TimePointsAvailable(i);
        tp_all = tp*Num_of_attributes;
        
        //X_tmp = X.row(i);
        //Z_tmp = Z.row(i);
        
        X_tmp = X(span(i), span(0, tp_all-1), span::all); //X.row(i);
        Z_tmp = Z(span(i), span(0, tp_all-1), span::all); //Z.row(i);

        if(updatealpha)
            Ri1_tmp = Ri_Within(i, TimePointsAvailable(i), alpha_samples.col(iter));
        else
            Ri1_tmp.eye(TimePointsAvailable(i), TimePointsAvailable(i));
        
        if(updatedelta)
            Ri0_tmp = Ri_Cross(i, Num_of_attributes, delta_samples.col(iter));
        else
            Ri0_tmp.eye(Num_of_attributes, Num_of_attributes);
        
        Ri_all_inv = kron(Ri1_tmp.i(), Ri0_tmp.i());
                
        Res = Y(span(0, tp_all-1), i) - Z_tmp* b_samples.slice(iter).row(i).t();
        
        //Rcout << "Res = " << Res << endl;
        
        //Rcout << "check 1" << endl;
        //Rcout << "X_tmp.t() * Ri_all_inv * X_tmp = " << X_tmp.t() * Ri_all_inv * X_tmp << endl;
        //Rcout << "Ind_mat_beta = " << Ind_mat_beta << endl;
        beta_cov += X_tmp.t() * Ri_all_inv * X_tmp + Ind_mat_beta/sigma2_beta;
        //Rcout << "check 2" << endl;
        //Rcout << "X_tmp.t() * Ri_all_inv * Res = " << X_tmp.t() * Ri_all_inv * Res;
        beta_mu  += X_tmp.t() * Ri_all_inv * Res;
    }
    //Rcout << "beta_cov = " << beta_cov << endl;
    //Rcout << "beta_mu = " << beta_mu << endl;
    beta_cov = inv_sympd(beta_cov);
    beta_mu = beta_cov*beta_mu;
    //Rcout << "beta_cov = " << beta_cov << endl;
    //Rcout << "beta_mu = " << beta_mu << endl;
    beta_samples.col(iter+1) = mvnrnd(beta_mu, beta_cov);
    
    //Rcout<< "Done Generate beta" <<endl;
}

void MultiLinearModel::Update_b(int iter)
{
    //Rcout << "Generate bi" << endl;
        
    mat Ri1_tmp, Ri0_tmp, Ri_all_inv;
    
    int bi_cov_size = Num_of_attributes*Num_of_RandomEffs;
    mat bi_cov = zeros<mat>(bi_cov_size, bi_cov_size);
    vec bi_mu = zeros<vec>(bi_cov_size);
    mat Z_tmp, X_tmp;
    vec Res;
    int tp, tp_all;
    
    for(int i=0; i<Num_of_obs; i++){
        tp = TimePointsAvailable(i);
        tp_all = tp*Num_of_attributes;
        
        if(updatealpha)
            Ri1_tmp = Ri_Within(i, TimePointsAvailable(i), alpha_samples.col(iter));
        else
            Ri1_tmp.eye(TimePointsAvailable(i), TimePointsAvailable(i));
        
        if(updatedelta)
            Ri0_tmp = Ri_Cross(i, Num_of_attributes, delta_samples.col(iter));
        else
            Ri0_tmp.eye(Num_of_attributes, Num_of_attributes);
        
        //Ri1_tmp = Ri(TimePointsAvailable(i), W(i), alpha_samples.col(iter));
        //Ri0_tmp = Ri(Num_of_attributes, V(i), delta_samples.col(iter));
        Ri_all_inv = kron(Ri1_tmp.i(), Ri0_tmp.i());
        //X_tmp = X.row(i);
        //Z_tmp = Z.row(i);
        X_tmp = X(span(i), span(0, tp_all-1), span::all); //X.row(i);
        Z_tmp = Z(span(i), span(0, tp_all-1), span::all); //Z.row(i);

        Res = Y(span(0, tp_all-1), i) - X_tmp* beta_samples.col(iter);
        bi_cov = inv_sympd( Z_tmp.t()*Ri_all_inv*Z_tmp + inv_sympd(Sigma_samples.slice(iter)) );
        bi_mu = bi_cov*Z_tmp.t()*Ri_all_inv*Res;
        b_samples.slice(iter+1).row(i) = mvnrnd(bi_mu, bi_cov).t();
    }
    //Rcout << "Done Generate bi" << endl;
}



void MultiLinearModel::Update_Sigma(int iter)
{
    //Rcout << "Generate Sigma" << endl;
    int bi_cov_size = Num_of_attributes*Num_of_RandomEffs;
    mat bi_cov = zeros<mat>(bi_cov_size, bi_cov_size);
    for(int i=0; i<Num_of_obs; i++)
        bi_cov += b_samples.slice(iter+1).row(i).t()*b_samples.slice(iter+1).row(i);
    bi_cov = bi_cov + Lambda;
    Sigma_samples.slice(iter+1) = iwishrnd( (bi_cov), (Num_of_obs + Vb));
    //Rcout << "End Sigma" << endl;
}

void MultiLinearModel::Update_alpha(int iter)
{
    //Rcout << "Generate alpha" << endl;
    double alpha_num=0, alpha_den=0;
    vec alpha_new;
    mat Ialpha_diag = eye(Num_of_alphas, Num_of_alphas);
    alpha_new =  mvnrnd(alpha_samples.col(iter), tuning_alpha*Ialpha_diag);
    
    //Rcout << "alpha_new = " << alpha_new << endl;
    
    mat Ri1_old, Ri1_new, Ri0;
    
    //Rcout << "Ri0 Done" << endl;
    
    mat Ri_inv_old, Ri_inv_new;
    
    mat X_tmp, Z_tmp;
    vec Res;
    int tp, tp_all;
    
    for(int i=0; i<Num_of_obs; i++){
        tp = TimePointsAvailable(i);
        tp_all = tp*Num_of_attributes;

        X_tmp = X(span(i), span(0, tp_all-1), span::all); //X.row(i);
        Z_tmp = Z(span(i), span(0, tp_all-1), span::all); //Z.row(i);

        
        Ri1_old = Ri_Within(i, TimePointsAvailable(i), alpha_samples.col(iter));
        
        Ri1_new = Ri_Within(i, TimePointsAvailable(i), alpha_new);
        
        //if(updatedelta)
             Ri0 = Ri_Cross(i, Num_of_attributes, delta_samples.col(iter));
        //else
         //   Ri0.eye(Num_of_attributes, Num_of_attributes);
        //Ri0 = Ri(Num_of_attributes, V(i), delta_samples.col(iter));
        
        
        if(0){
            Rcout << "W(i) = \n" << W(i).row(0) << endl;
            Rcout << "W(i) = \n" << W(i).row(1) << endl;
            Rcout << "alpha_new = \n" << alpha_new << endl;
            Rcout << " Ri1_new = \n" << Ri1_new << endl;
        }
        
        Ri_inv_old = kron(Ri1_old.i(), Ri0.i());
        Ri_inv_new = kron(Ri1_new.i(), Ri0.i());
        
        //Rcout << "size (X_tmp) = " << size(X_tmp) << endl;
        //Rcout << "size of b = " << size(b_samples.slice(iter+1).row(i)) << endl;
        
        Res = Y(span(0, tp_all-1), i) - X_tmp* beta_samples.col(iter+1)-Z_tmp* b_samples.slice(iter+1).row(i).t();
 
        alpha_num += -0.5 * Num_of_attributes* log(det(Ri1_new)) - 0.5*as_scalar(Res.t()*Ri_inv_new*Res);
        alpha_den += -0.5 * Num_of_attributes* log(det(Ri1_old)) - 0.5*as_scalar(Res.t()*Ri_inv_old*Res);
        
    }
    
    
    //vec bb = { -1.2, -0.5};
    //Rcout << "Ri1 = \n" << Ri_Within(0, 4, bb);

    
    alpha_num += -0.5*accu(square(alpha_new))/sigma2_alpha;
    alpha_den += -0.5*accu(square(alpha_samples.col(iter)))/sigma2_alpha;
    
    alpha_samples.col(iter+1) = alpha_samples.col(iter);
    if(log(Rf_runif(0., 1.)) < alpha_num - alpha_den ){
        alpha_samples.col(iter+1) = alpha_new;
        acc_rate_alpha++;
    }
    
    if((iter+1)%500 == 0){
        if( acc_rate_alpha/iter<0.25 )
            tuning_alpha = tuning_alpha/2.;
        if( (1.*acc_rate_alpha)/iter>0.50 )
            tuning_alpha = 2*tuning_alpha;
        
    }
     
    //Rcout << "Done Generate alpha" << endl;
}

void MultiLinearModel::Update_delta(int iter)
{
    //Rcout << "Generate delta" << endl;
    double delta_num=0, delta_den=0;
    vec delta_new;
    mat Idelta_diag = eye(Num_of_deltas, Num_of_deltas);
    delta_new =  mvnrnd(delta_samples.col(iter), tuning_delta*Idelta_diag);
    
    //delta_new = delta_samples.col(iter);
    //delta_samples.col(iter+1) = delta_samples.col(iter);
    
    mat Ri0_old, Ri0_new;
        
    mat Ri1;
    
    mat Ri_inv_old, Ri_inv_new;
    
    mat X_tmp, Z_tmp;
    vec Res;
    int tp, tp_all;
    
    for(int i=0; i<Num_of_obs; i++){
        //delta_new(j) = delta_samples(j, iter+1) + tuning_delta*randn(); //randn( distr_param(delta_samples(j, i),tuning_delta) );
        tp = TimePointsAvailable(i);
        tp_all = tp*Num_of_attributes;
        X_tmp = X(span(i), span(0, tp_all-1), span::all); //X.row(i);
        Z_tmp = Z(span(i), span(0, tp_all-1), span::all); //Z.row(i);

        Ri1 = Ri_Within(i, TimePointsAvailable(i), alpha_samples.col(iter+1));
        //if(updatealpha)
        //    Ri1 = Ri_Within(i, TimePointsAvailable(i), alpha_samples.col(iter+1));
        //else
        //    Ri1.eye(TimePointsAvailable(i), TimePointsAvailable(i));
        
        
        //Ri1 = Ri(TimePointsAvailable(i), W(i), alpha_samples.col(iter+1));
        Ri0_old = Ri_Cross(i, Num_of_attributes, delta_samples.col(iter));
        Ri0_new = Ri_Cross(i, Num_of_attributes, delta_new);
        
        if(0){
            Rcout << "V(i) = \n" << V(i).row(0) << endl;
            Rcout << "V(i) = \n" << V(i).row(1) << endl;
            Rcout << "V(i) = \n" << V(i).row(2) << endl;
            Rcout << "delta_new = \n" << delta_new << endl;
            Rcout << " Ri0_new = \n" << Ri0_new << endl;
        }
        
        Ri_inv_old = kron(Ri1.i(), Ri0_old.i());
        Ri_inv_new = kron(Ri1.i(), Ri0_new.i());
        
        Res = Y(span(0, tp_all-1), i) - X_tmp* beta_samples.col(iter+1)-Z_tmp* b_samples.slice(iter+1).row(i).t();
 
        delta_num += -0.5 * TimePointsAvailable(i)* log(det(Ri0_new)) - 0.5*as_scalar(Res.t()*Ri_inv_new*Res);
        delta_den += -0.5 * TimePointsAvailable(i)* log(det(Ri0_old)) - 0.5*as_scalar(Res.t()*Ri_inv_old*Res);
        
    }
    
    //Rcout << "delta_num = " << delta_num << endl;
    //Rcout << "delta_den = " << delta_den << endl;

    //Rcout << " delta_new = " << endl << delta_new << endl;
    //Rcout << " delta_old = " << endl << delta_samples.col(iter) << endl;
    
    
    
    delta_num += -0.5*accu(square(delta_new))/sigma2_delta;
    delta_den += -0.5*accu(square(delta_samples.col(iter)))/sigma2_delta;
    
    //Rcout << "delta_num = " << delta_num << endl;
    //Rcout << "delta_den = " << delta_den << endl;

    delta_samples.col(iter+1) = delta_samples.col(iter);
    if(log(Rf_runif(0., 1.)) < delta_num - delta_den){
        //Rcout << "j = " << j << endl;
        delta_samples.col(iter+1) = delta_new;
        //if(j==0)
            acc_rate_delta++;
    }
    //else delta_new(j) = delta_samples(j, iter+1);
    
    //Rcout << "acc_rate_delta = " << acc_rate_delta << endl;
    //vec aa = {-0.5, -0.4, -0.8, -0.2, -0.3};
    //Rcout << "Ri0 = \n" << Ri_Cross(0, 5, aa);
    
    //Rcout << "===============================" << endl;
    if((iter+1)%500 == 0){
        if( acc_rate_delta/iter<0.25 )
            tuning_delta = tuning_delta/2.;
        if( (1.*acc_rate_delta)/iter>0.50 )
            tuning_delta = 2*tuning_delta;
        
    }
        
    
    //Rcout << "Done Generate delta" << endl;
    
}



 void MultiLinearModel::Update_delta_CompWise(int iter)
 {
     //Rcout << "Generate Update_delta_CompWise" << endl;
     double delta_num=0, delta_den=0;
     vec delta_new;
     mat Idelta_diag = eye(Num_of_deltas, Num_of_deltas);
     //delta_new =  mvnrnd(delta_samples.col(iter), tuning_delta*Idelta_diag);
     
     delta_new = delta_samples.col(iter);
     delta_samples.col(iter+1) = delta_samples.col(iter);
     
     mat Ri0_old, Ri0_new;
         
     mat Ri1;
     
     mat Ri_inv_old, Ri_inv_new;
     
     mat X_tmp, Z_tmp;
     vec Res;
     int tp, tp_all;
     
     for(int  j=0; j<Num_of_deltas; j++){
         
         for(int i=0; i<Num_of_obs; i++){
             delta_new(j) = delta_samples(j, iter+1) + tuning_delta*randn(); //randn( distr_param(delta_samples(j, i),tuning_delta) );
             tp = TimePointsAvailable(i);
             tp_all = tp*Num_of_attributes;
             X_tmp = X(span(i), span(0, tp_all-1), span::all); //X.row(i);
             Z_tmp = Z(span(i), span(0, tp_all-1), span::all); //Z.row(i);

             if(updatealpha)
                 Ri1 = Ri_Within(i, TimePointsAvailable(i), alpha_samples.col(iter+1));
             else
                 Ri1.eye(TimePointsAvailable(i), TimePointsAvailable(i));
             
             
             //Ri1 = Ri(TimePointsAvailable(i), W(i), alpha_samples.col(iter+1));
             Ri0_old = Ri_Cross(i, Num_of_attributes, delta_samples.col(iter));
             Ri0_new = Ri_Cross(i, Num_of_attributes, delta_new);
             
             if(0){
                 Rcout << "V(i) = \n" << V(i).row(0) << endl;
                 Rcout << "V(i) = \n" << V(i).row(1) << endl;
                 Rcout << "V(i) = \n" << V(i).row(2) << endl;
                 Rcout << "delta_new = \n" << delta_new << endl;
                 Rcout << " Ri0_new = \n" << Ri0_new << endl;
             }
             
             Ri_inv_old = kron(Ri1.i(), Ri0_old.i());
             Ri_inv_new = kron(Ri1.i(), Ri0_new.i());
             
             Res = Y(span(0, tp_all-1), i) - X_tmp* beta_samples.col(iter+1)-Z_tmp* b_samples.slice(iter+1).row(i).t();
      
             delta_num += -0.5 * TimePointsAvailable(i)* log(det(Ri0_new)) - 0.5*as_scalar(Res.t()*Ri_inv_new*Res);
             delta_den += -0.5 * TimePointsAvailable(i)* log(det(Ri0_old)) - 0.5*as_scalar(Res.t()*Ri_inv_old*Res);
             
         }
         
         //Rcout << "delta_num = " << delta_num << endl;
         //Rcout << "delta_den = " << delta_den << endl;

         //Rcout << " delta_new = " << endl << delta_new << endl;
         //Rcout << " delta_old = " << endl << delta_samples.col(iter) << endl;
         
         
         
         delta_num += -0.5*accu(square(delta_new))/sigma2_delta;
         delta_den += -0.5*accu(square(delta_samples.col(iter)))/sigma2_delta;
         
         //Rcout << "delta_num = " << delta_num << endl;
         //Rcout << "delta_den = " << delta_den << endl;

         if(log(Rf_runif(0., 1.)) < delta_num - delta_den){
             //Rcout << "j = " << j << endl;
             delta_samples(j, iter+1) = delta_new(j);
             if(j==0)
                 acc_rate_delta++;
         }
         else delta_new(j) = delta_samples(j, iter+1);
         
         if((iter+1)%100==0){
             Rcout << "acc_rate_delta = " << acc_rate_delta << endl;
         //vec aa = {-0.5, -0.4, -0.8};
         //Rcout << "Ri0 = \n" << Ri_Cross(0, 3, aa);
         
             Rcout << "===============================" << endl;
         }
         if((iter+1)%500 == 0){
             if( acc_rate_delta/iter<0.25 )
                 tuning_delta = tuning_delta/2.;
             if( (1.*acc_rate_delta)/iter>0.50 )
                 tuning_delta = 2*tuning_delta;
             
         }
         
     }
     //Rcout << "Done Generate delta" << endl;
     
 }


void MultiLinearModel::ParameterEstimation()
{
    //Rcout << "Time = " << TimePointsAvailable << endl;
    beta_mean = mean(beta_samples, 1);
    if(updatealpha)
        alpha_mean = mean(alpha_samples, 1);
    if(updatedelta)
       delta_mean = mean(delta_samples, 1);
    Sigma_mean = mean(Sigma_samples, 2);
    b_mean = mean(b_samples, 2);

    //Rcout << "beta est = " << endl << beta_mean << endl;
    //Rcout << "alpha est = " << endl << alpha_mean << endl;
    //Rcout << "delta est = " << endl << delta_mean << endl;
    //Rcout << "Sigma est = " << endl << Sigma_mean << endl;
    
    mat X_tmp, Z_tmp;
    vec y_mu;
    int tp, tp_all;
    vec pit;
    uvec diff_index;
    logL = 0.;
    AIC = 0.;
    BIC = 0.;
    DIC = 0.;
    MPL = 0.;
    ePCP = 0.; 
    vec CPO_tmp;
    mat CPO = zeros<mat>(Num_of_obs, Num_of_attributes*TimePointsAvailable.max());
    mat Rho11, Rho21, Rho22, Rho1, Rho2, Rho1_inv;
    double rho_cAIC = 0.;
    
    double PMC =0.;
    uvec num_of_ones;

    
    for(int i=0; i<Num_of_obs; i++){

        tp = TimePointsAvailable(i);
        tp_all = tp*Num_of_attributes;


        X_tmp = X(span(i), span(0, tp_all-1), span::all); //X.row(i);
        Z_tmp = Z(span(i), span(0, tp_all-1), span::all); //Z.row(i);
                
        y_mu = X_tmp* beta_mean+Z_tmp*b_mean.row(i).t();
        
        //Rcout << "Y_Binary(span(0, tp_all-1), i) = " << Y_Binary(span(0, tp_all-1), i) << endl;
        
        PCP += 1.*accu(1.*(y_mu>0) == Y_Binary(span(0, tp_all-1), i));
        
        
        pit = normcdf(y_mu, 0., 1.);
        //Rcout << "i=" << i << ", t=" << t << "\tpit=" << pit << endl;
        
        diff_index = find( (pit-Y_Binary(span(0, tp_all-1), i))!=0.);
        
        
        logL += accu(Y_Binary(span(0, tp_all-1), i)%log(pit) + (1-Y_Binary(span(0, tp_all-1), i))%log(1-pit));
        
        ePCP += accu(Y_Binary(span(0, tp_all-1), i)%pit + (1-Y_Binary(span(0, tp_all-1), i))%(1-pit));
        
        //Rcout << "diff_index = " << diff_index << endl;

        
        for(int iter = Num_of_iterations/2; iter<Num_of_iterations; iter++){
            pit = normcdf( (X_tmp*beta_samples.col(iter) + Z_tmp*b_samples.slice(iter).row(i).t()) , 0., 1.);
            
            diff_index = find( (pit-Y_Binary(span(0, tp_all-1), i))!=0.);
            
            CPO_tmp = (Y_Binary(span(0, tp_all-1), i)%log(pit) + (1-Y_Binary(span(0, tp_all-1), i))%log(1-pit));
            
            DIC += accu(CPO_tmp);
  
            //Rcout << "exp(-CPO_tmp) = " << exp(-CPO_tmp) << endl;
            CPO(i, span(0, tp_all-1)) += exp(-CPO_tmp).t();
        }

        //Rcout << "Check 11" << endl;
        //Rho11 = X_tmp.submat(0, 0, tp-1, Num_of_FixedEffs-1).t() * X_tmp.submat(0, 0, tp-1, Num_of_FixedEffs-1);
        //Rho21 = Z_tmp.submat(0, 0, tp-1, Num_of_RandomEffs-1).t() * X_tmp.submat(0, 0, tp-1, Num_of_FixedEffs-1);
        //Rho22 = Z_tmp.submat(0, 0, tp-1, Num_of_RandomEffs-1).t() * Z_tmp.submat(0, 0, tp-1, Num_of_RandomEffs-1);
        
        Rho11 = X_tmp.t() * X_tmp;
        Rho21 = Z_tmp.t() * X_tmp;
        Rho22 = Z_tmp.t() * Z_tmp;

        
        //Rcout << "Rho11 = " << Rho11 << endl;
        //Rcout << "Rho21 = " << Rho21 << endl;
        //Rcout << "Rho22 = " << Rho22 << endl;
        
        //Rcout << "Check 12" << endl;
        
        //Rcout << "join_rows(Rho11, Rho12) = " << endl <<  join_cols(Rho11, Rho21) << endl;
        //Rcout << "join_rows(Rho12.t(), Rho22) = " << endl <<  join_rows(Rho12.t(), Rho22) << endl;
        
        Rho1 = join_rows(join_cols(Rho11, Rho21), join_cols(Rho21.t(), Rho22+Sigma_mean)); //
        
        //Rcout << "size(Rho1)=" << size(Rho1) << endl;
        
        //Rcout << "det( Rho11 ) =" << det(Rho11) << endl;
        
        if(det(Rho1)==0){
            Rho1_inv = Rho1+ 0.01*eye(size(Rho1));
            Rho1_inv = Rho1_inv.i();
        }
        else
            Rho1_inv = Rho1.i();
        //Rcout << "size(Rho11)=" << size(Rho11) << endl;
        //Rcout << Rho11 << endl;
        //Rcout <<  Rho11.i() << endl;
        
        Rho2 = join_rows(join_cols(Rho11, Rho21), join_cols(Rho21.t(), Rho22));
        //Rcout << "Check 14" << endl;
        //Rcout << "Rho1 = " << Rho1 << endl;
        //Rcout << "Rho1.i() = " << Rho11_adj.i() << endl;
        rho_cAIC += trace(Rho1_inv * Rho2);
        
        //Rcout << "Check 15" << endl;
            
    }

    
    
    
    CPO = 1./CPO;
    
    //Rcout << "CPO = " << endl << CPO.submat(0, 0, 9, 3) << endl;
    
    CPO.elem( find_nonfinite(CPO) ).zeros();
    
    MPL = accu(CPO);
    
    PCP = PCP/(accu(TimePointsAvailable)*Num_of_attributes);
    ePCP = ePCP/(accu(TimePointsAvailable)*Num_of_attributes);
    
    
    num_of_ones = find(Y_Binary==1);
    PMC = num_of_ones.n_elem/(accu(TimePointsAvailable)*Num_of_attributes);
    
    //Rcout << "Num.of.ones = " << num_of_ones.n_elem << endl;

    
    //Rcout << "PMC = " << PMC << endl;
    
    PRE = (PCP-PMC) / (1.-PMC);
        
    AIC = -2*logL + 2 * Num_of_attributes*(Num_of_FixedEffs+Num_of_obs*Num_of_RandomEffs) + (Num_of_deltas+Num_of_alphas);
    
    cAIC = -2*logL +2 * (rho_cAIC+1);
    
    BIC = -2*logL + log(Num_of_obs) * ( Num_of_attributes*(Num_of_FixedEffs+Num_of_obs*Num_of_RandomEffs)+ (Num_of_deltas+Num_of_alphas));
    
    DIC = -4*DIC/(Num_of_iterations/2) + 2*logL;
    
    Rcout << "PCP = " << PCP << endl;
    Rcout << "PRE = " << PRE << endl;
    Rcout << "ePCP = " << ePCP << endl;
    Rcout << "logL = " << logL << endl;
    Rcout << "AIC = " << AIC << endl;
    Rcout << "cAIC = " << cAIC << endl;
    Rcout << "BIC = " << BIC << endl;
    Rcout << "DIC = " << DIC << endl;
    Rcout << "MPL = " << MPL << endl;

    
/*
    rowvec X_tmp, Z_tmp, Ri_tmp;

    double pit, CPO_tmp, ESS_GP_tmp, RJ1, RJ2, ESS=0., GP=0.;
    logL = 0.;
    mat Djt(Num_of_Timepoints, Num_of_covariates, fill::zeros);
    mat Omega_I(Num_of_covariates, Num_of_covariates, fill::zeros), M_LZ(Num_of_covariates, Num_of_covariates, fill::zeros);
    vec mu_it(Num_of_Timepoints, fill::zeros), p_it(Num_of_Timepoints, fill::zeros);
    mat A_sqrt, Cov_Y, V, V_inv, Omega_I_inv, V_LZ, Gamma_RJ, Djt_sub;
    
    mat CPO = zeros<mat>(Num_of_obs, TimePointsAvailable.max());
    int tp;
    
    CIC = 0.;
    RJ_R = 0.;
    ACC = 0.;

    //Rcout << "============== 0 ============"<<endl;
    for(int i=0; i<Num_of_obs; i++){
        //Rcout << "i = " << i << endl;
        tp = TimePointsAvailable(i);
        for(int t=0; t<tp; t++){
            
            X_tmp = X.slice(i).row(t);
            Z_tmp = Z.slice(i).row(t);
            //Rcout << "X_tmp = " << X_tmp << endl;
            //Rcout << "Z_tmp = " << Z_tmp << endl;
            mu_it(t) = as_scalar( X_tmp*beta_mean + Z_tmp*b_mean.col(i));
            p_it(t) = normcdf(mu_it(t));
            
            ACC += 1.*(1.*(mu_it(t)>0) == Y(t, i));
            
            for(int j=0; j<Num_of_covariates; j++)
                Djt(t, j) = X(t,j,i)/normpdf(Rf_pnorm5(mu_it(t), 0., 1., 1, 0));
        }
        //Rcout << "============== 1 ============"<<endl;
        Djt_sub = Djt.head_rows(tp);
        //Rcout << "============== 2 ============"<<endl;
        Cov_Y = (Y(span(0, tp-1), i)-p_it.head(tp))*(Y(span(0, tp-1), i)-p_it.head(tp)).t();
        //Rcout << "============== 3 ============" << endl;
        A_sqrt = diagmat(sqrt(p_it.head(tp)));
        if(updatedelta)
            V = A_sqrt*Ri_Version2(i, tp, delta_mean)*A_sqrt;
        else
            V = A_sqrt*A_sqrt;
        //Rcout << "============== 4 ============"<<endl;
        V_inv = V.i();
        Omega_I += Djt_sub.t()*V_inv*Djt_sub;
        //Rcout << "============== 5 ============"<<endl;
        M_LZ += Djt_sub.t()*V_inv*Cov_Y*V_inv*Djt_sub;
        //Rcout << "============== 6 ============"<<endl;
        ESS_GP_tmp = as_scalar((Y(span(0, tp-1), i)-mu_it.head(tp)).t()*V_inv*(Y(span(0, tp-1), i)-mu_it.head(tp)));
        //Rcout << "============== 7 ============"<<endl;
        ESS += ESS_GP_tmp;
        GP += -0.5*ESS_GP_tmp + log(det(V));
    }
    //Rcout << "============== 2 ============"<<endl;
    //Rcout << " Omega_I = " << endl << Omega_I << endl;
    Omega_I_inv = Omega_I.i();
    //Rcout << "============== 3 ============"<<endl;
    Gamma_RJ = Omega_I_inv*M_LZ;
    //Rcout << "============== 4 ============"<<endl;
    V_LZ = Gamma_RJ*Omega_I_inv;
    //Rcout << "============== 5 ============"<<endl;
    RJ1 = trace(Gamma_RJ)/Num_of_covariates;
    RJ2 = accu((diagvec(Gamma_RJ*Gamma_RJ)))/Num_of_covariates;

    RJ_R = sqrt((1-RJ1)*(1-RJ1)+(1-RJ2)*(1-RJ2));
    //SC = ESS/(N-P-a
    //GP = -0.5*GP

    CIC = trace(Omega_I*V_LZ);
    //Rcout << "RJ_R = " << RJ_R << "\tCIC = " << CIC << endl;
    //cat("RJ.R = ", RJ.R, "\t", "SC = ", SC, "\n")

    
    for(int i=0; i<Num_of_obs; i++){
        for(int t=0; t<TimePointsAvailable(i); t++){
            
            X_tmp = X.slice(i).row(t);
            Z_tmp = Z.slice(i).row(t);
                    
            
            for(int iter = Num_of_iterations/2; iter<Num_of_iterations; iter++){
                pit = normcdf(as_scalar( X_tmp*beta_samples.col(iter) + Z_tmp*b_samples.slice(iter).col(i) ) , 0., 1.);
                
                if(pit == 1 && Y(t, i) == 1){
                    DIC += 0.;
                    CPO_tmp = 0.;
                }
                else if(pit == 0 && Y(t, i) == 0){
                    DIC += 0.;
                //else if(pit == 0 && Y(t, i) == 1)
                //    Likelihood += 0.;
                //else if(pit == 1 && Y(t, i) == 0)
                //    Likelihood += 0.;
                    CPO_tmp = 0.;
                }
                else{
                    //Likelihood *= pow(pit, Y(t, i))*pow( (1-pit), (1-Y(t, i)) );
                    CPO_tmp = Y(t, i)*log(pit) + (1-Y(t, i))*log(1-pit);
                    DIC += CPO_tmp;
                    
                }
                
                CPO(i, t) += exp(-CPO_tmp);
            }
    
            
            
            pit = normcdf(as_scalar( X_tmp*beta_mean + Z_tmp*b_mean.col(i) ), 0., 1.);
            //Rcout << "i=" << i << ", t=" << t << "\tpit=" << pit << endl;
            
            if(pit == 1 && Y(t, i) == 1)
                logL += 0.;
            else if(pit == 0 && Y(t, i) == 0)
                logL += 0.;
            //else if(pit == 0 && Y(t, i) == 1)
            //    Likelihood += 0.;
            //else if(pit == 1 && Y(t, i) == 0)
            //    Likelihood += 0.;
            else
                //Likelihood *= pow(pit, Y(t, i))*pow( (1-pit), (1-Y(t, i)) );
                logL += Y(t, i)*log(pit) + (1-Y(t, i))*log(1-pit);
        }
    }

    CPO = 1./CPO;
    
    //Rcout << "CPO = " << endl << CPO.submat(0, 0, 9, 3) << endl;
    
    CPO.elem( find_nonfinite(CPO) ).zeros();
    
    MPL = accu(CPO);
    DIC = -4*DIC/(Num_of_iterations/2) + 2*logL;
    AIC = -2*logL + 2 * (Num_of_covariates+Num_of_obs*Num_of_RanEffs + Num_of_deltas);
    BIC = -2*logL + log(Num_of_obs) * (Num_of_covariates+Num_of_obs*Num_of_RanEffs+ Num_of_deltas);
*/
}



SEXP MultiLinearModel::MCMC_Procedure()
{
    Rcout << "============= FMR: MCMC Starts=============="<< endl;
    List PosteriorSamples;
    List PosteriorEstimates;
    List Posterior;
    
    time_t start = time(NULL);
    
    int iter = 0;
    
    while(iter < Num_of_iterations-1){
        if(updateystar)
            Update_ystar(iter);
        
        if(updatebeta)
            Update_beta(iter);
        else
            beta_samples.col(iter+1) = beta_samples.col(iter);
        
        if(updateb)
            Update_b(iter);
        else
            b_samples.slice(iter+1) = b_samples.slice(iter);
        
        if(updateSigma)
            Update_Sigma(iter);
        else
            Sigma_samples.slice(iter+1) = Sigma_samples.slice(iter);
        
        //Rcout << "Check 1 update alpha" << updatealpha << endl;
       
        if(updatealpha)
            Update_alpha(iter);
        //else
            //alpha_samples.col(iter+1) = alpha_samples.col(iter);
        
        //Rcout << "Check 2 update delta" << updatealpha << endl;
        
        if(updatedelta){
            if(Num_of_deltas>3)
                Update_delta_CompWise(iter);
            else
                Update_delta(iter);
        }
        //else
            //delta_samples.col(iter+1) = delta_samples.col(iter);

        iter++;
        if((iter+1)%100==0)
            Rcout << iter+1 << endl;
        //Rcout<<"\t"<< round((iter+1.)/Num_of_iterations*100)<<" %"<<'\r';
    }
    ParameterEstimation();
        
    Rcout << "============= FMR: MCMC: Done =============="<< endl;
    time_t end = time(NULL);
    Rcout<<"Execution Time: "<< (double)(end-start)<<" Seconds"<<std::endl;
    
    PosteriorSamples["beta.samples"] = beta_samples;
    PosteriorSamples["b.samples"] = b_samples;
    PosteriorSamples["Sigma.samples"] = Sigma_samples;
    if(updatealpha)
        PosteriorSamples["alpha.samples"] = alpha_samples;
    if(updatedelta)
        PosteriorSamples["delta.samples"] = delta_samples;

    PosteriorEstimates["beta.mean"] = beta_mean;
    if(updatealpha)
        PosteriorEstimates["alpha.mean"] = alpha_mean;
    if(updatedelta)
        PosteriorEstimates["delta.mean"] = delta_mean;
    PosteriorEstimates["b.mean"] = b_mean;
    PosteriorEstimates["Sigma.mean"] = Sigma_mean;

    
    PosteriorEstimates["AIC"] = AIC;
    PosteriorEstimates["cAIC"] = cAIC;
    PosteriorEstimates["BIC"] = BIC;
    PosteriorEstimates["logL"] = logL;
    PosteriorEstimates["DIC"] = DIC;
    PosteriorEstimates["MPL"] = MPL;
    PosteriorEstimates["PCP"] = PCP;
    PosteriorEstimates["PRE"] = PRE;
    PosteriorEstimates["ePCP"] = ePCP;
    PosteriorEstimates["Ri0"] = Ri_Cross(0, Num_of_attributes, delta_mean);
    PosteriorEstimates["Ri1"] = Ri_Within(0, TimePointsAvailable.max(), alpha_mean);
    
    Posterior["PosteriorEstimates"] = PosteriorEstimates;
    Posterior["PosteriorSamples"] = PosteriorSamples;
   
    return (Posterior);
}


