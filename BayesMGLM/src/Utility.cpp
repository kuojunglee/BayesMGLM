#include "MLModelSelection.h"
#include "MGLM.h"

// [[Rcpp::export]]
RcppExport SEXP MLModelSelectionMCMC(SEXP i_Num_of_iterations, SEXP list_Data, SEXP list_InitialValues, SEXP list_HyperPara, SEXP list_UpdatePara, SEXP list_TuningPara)
{
    List lData(list_Data);
    List lInitialValues(list_InitialValues);
    List lHyperPara(list_HyperPara);
    List lUpdatePara(list_UpdatePara);
    List lTuningPara(list_TuningPara);
    
    List PosteriorSamples;
    
    int iNum_of_iterations = Rcpp::as<int> (i_Num_of_iterations);
    
    MLModelSelection DoMLModelSelectionMCMC(iNum_of_iterations, lData, lInitialValues, lHyperPara, lUpdatePara, lTuningPara);
    
    PosteriorSamples = DoMLModelSelectionMCMC.MCMC_Procedure();
    
    //Rcout << "Check1" << endl;
    return PosteriorSamples;
    
}

// [[Rcpp::export]]
RcppExport SEXP MultiLinearModelMCMC(SEXP i_Num_of_iterations, SEXP list_Data, SEXP list_InitialValues, SEXP list_HyperPara, SEXP list_UpdatePara, SEXP list_TuningPara)
{
    List lData(list_Data);
    List lInitialValues(list_InitialValues);
    List lHyperPara(list_HyperPara);
    List lUpdatePara(list_UpdatePara);
    List lTuningPara(list_TuningPara);

    List PosteriorSamples;
    
    int iNum_of_iterations = Rcpp::as<int> (i_Num_of_iterations);
    
    MultiLinearModel DoMultiLinearModelMCMC(iNum_of_iterations, lData, lInitialValues, lHyperPara, lUpdatePara, lTuningPara);
    
    PosteriorSamples = DoMultiLinearModelMCMC.MCMC_Procedure();
    
    //Rcout << "Check1" << endl;
    return PosteriorSamples;
    
}

