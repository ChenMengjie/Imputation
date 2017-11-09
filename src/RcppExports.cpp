// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// calculate_weights
arma::vec calculate_weights(arma::vec z, arma::mat X);
RcppExport SEXP _Imputation_calculate_weights(SEXP zSEXP, SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(calculate_weights(z, X));
    return rcpp_result_gen;
END_RCPP
}
// fitting_lasso
Rcpp::List fitting_lasso(arma::vec y, arma::mat X, bool min);
RcppExport SEXP _Imputation_fitting_lasso(SEXP ySEXP, SEXP XSEXP, SEXP minSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< bool >::type min(minSEXP);
    rcpp_result_gen = Rcpp::wrap(fitting_lasso(y, X, min));
    return rcpp_result_gen;
END_RCPP
}
// log_factorial
double log_factorial(int Y);
RcppExport SEXP _Imputation_log_factorial(SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(log_factorial(Y));
    return rcpp_result_gen;
END_RCPP
}
// log_factorial_calculated
arma::vec log_factorial_calculated(int N);
RcppExport SEXP _Imputation_log_factorial_calculated(SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(log_factorial_calculated(N));
    return rcpp_result_gen;
END_RCPP
}
// Mix_gradient_and_LogLikelihood_for_individual_sample
Rcpp::List Mix_gradient_and_LogLikelihood_for_individual_sample(arma::vec Y, arma::vec W, arma::vec V, arma::vec WY, arma::vec WWY, arma::vec W3Y, arma::vec W4Y, arma::vec VY, arma::vec VVY, arma::vec V3Y, arma::vec V4Y, arma::vec WW, arma::vec W3, arma::vec W4, arma::vec VV, arma::vec V3, arma::vec V4, double a0, double a1, double a2, double a3, double a4, double b1, double b2, double b3, double b4, double psi, int n, double sum_log_factorial_Y);
RcppExport SEXP _Imputation_Mix_gradient_and_LogLikelihood_for_individual_sample(SEXP YSEXP, SEXP WSEXP, SEXP VSEXP, SEXP WYSEXP, SEXP WWYSEXP, SEXP W3YSEXP, SEXP W4YSEXP, SEXP VYSEXP, SEXP VVYSEXP, SEXP V3YSEXP, SEXP V4YSEXP, SEXP WWSEXP, SEXP W3SEXP, SEXP W4SEXP, SEXP VVSEXP, SEXP V3SEXP, SEXP V4SEXP, SEXP a0SEXP, SEXP a1SEXP, SEXP a2SEXP, SEXP a3SEXP, SEXP a4SEXP, SEXP b1SEXP, SEXP b2SEXP, SEXP b3SEXP, SEXP b4SEXP, SEXP psiSEXP, SEXP nSEXP, SEXP sum_log_factorial_YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type V(VSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type WY(WYSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type WWY(WWYSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type W3Y(W3YSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type W4Y(W4YSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type VY(VYSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type VVY(VVYSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type V3Y(V3YSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type V4Y(V4YSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type WW(WWSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type W3(W3SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type W4(W4SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type VV(VVSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type V3(V3SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type V4(V4SEXP);
    Rcpp::traits::input_parameter< double >::type a0(a0SEXP);
    Rcpp::traits::input_parameter< double >::type a1(a1SEXP);
    Rcpp::traits::input_parameter< double >::type a2(a2SEXP);
    Rcpp::traits::input_parameter< double >::type a3(a3SEXP);
    Rcpp::traits::input_parameter< double >::type a4(a4SEXP);
    Rcpp::traits::input_parameter< double >::type b1(b1SEXP);
    Rcpp::traits::input_parameter< double >::type b2(b2SEXP);
    Rcpp::traits::input_parameter< double >::type b3(b3SEXP);
    Rcpp::traits::input_parameter< double >::type b4(b4SEXP);
    Rcpp::traits::input_parameter< double >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type sum_log_factorial_Y(sum_log_factorial_YSEXP);
    rcpp_result_gen = Rcpp::wrap(Mix_gradient_and_LogLikelihood_for_individual_sample(Y, W, V, WY, WWY, W3Y, W4Y, VY, VVY, V3Y, V4Y, WW, W3, W4, VV, V3, V4, a0, a1, a2, a3, a4, b1, b2, b3, b4, psi, n, sum_log_factorial_Y));
    return rcpp_result_gen;
END_RCPP
}
// Mix_LogLikelihood_for_individual_sample
double Mix_LogLikelihood_for_individual_sample(arma::vec Y, arma::vec W, arma::vec V, arma::vec WW, arma::vec VV, arma::vec W3, arma::vec V3, arma::vec W4, arma::vec V4, double a0, double a1, double a2, double a3, double a4, double b1, double b2, double b3, double b4, double psi, int n, double sum_log_factorial_Y);
RcppExport SEXP _Imputation_Mix_LogLikelihood_for_individual_sample(SEXP YSEXP, SEXP WSEXP, SEXP VSEXP, SEXP WWSEXP, SEXP VVSEXP, SEXP W3SEXP, SEXP V3SEXP, SEXP W4SEXP, SEXP V4SEXP, SEXP a0SEXP, SEXP a1SEXP, SEXP a2SEXP, SEXP a3SEXP, SEXP a4SEXP, SEXP b1SEXP, SEXP b2SEXP, SEXP b3SEXP, SEXP b4SEXP, SEXP psiSEXP, SEXP nSEXP, SEXP sum_log_factorial_YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type V(VSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type WW(WWSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type VV(VVSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type W3(W3SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type V3(V3SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type W4(W4SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type V4(V4SEXP);
    Rcpp::traits::input_parameter< double >::type a0(a0SEXP);
    Rcpp::traits::input_parameter< double >::type a1(a1SEXP);
    Rcpp::traits::input_parameter< double >::type a2(a2SEXP);
    Rcpp::traits::input_parameter< double >::type a3(a3SEXP);
    Rcpp::traits::input_parameter< double >::type a4(a4SEXP);
    Rcpp::traits::input_parameter< double >::type b1(b1SEXP);
    Rcpp::traits::input_parameter< double >::type b2(b2SEXP);
    Rcpp::traits::input_parameter< double >::type b3(b3SEXP);
    Rcpp::traits::input_parameter< double >::type b4(b4SEXP);
    Rcpp::traits::input_parameter< double >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type sum_log_factorial_Y(sum_log_factorial_YSEXP);
    rcpp_result_gen = Rcpp::wrap(Mix_LogLikelihood_for_individual_sample(Y, W, V, WW, VV, W3, V3, W4, V4, a0, a1, a2, a3, a4, b1, b2, b3, b4, psi, n, sum_log_factorial_Y));
    return rcpp_result_gen;
END_RCPP
}
// Mix_select_stepsize_for_a_parameter
double Mix_select_stepsize_for_a_parameter(arma::vec Y, arma::vec W, arma::vec V, arma::vec WW, arma::vec VV, arma::vec W3, arma::vec V3, arma::vec W4, arma::vec V4, double ll, double sum_log_factorial_Y, arma::vec gradient, arma::vec parameters, int ind, double gamma, int n, double down);
RcppExport SEXP _Imputation_Mix_select_stepsize_for_a_parameter(SEXP YSEXP, SEXP WSEXP, SEXP VSEXP, SEXP WWSEXP, SEXP VVSEXP, SEXP W3SEXP, SEXP V3SEXP, SEXP W4SEXP, SEXP V4SEXP, SEXP llSEXP, SEXP sum_log_factorial_YSEXP, SEXP gradientSEXP, SEXP parametersSEXP, SEXP indSEXP, SEXP gammaSEXP, SEXP nSEXP, SEXP downSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type V(VSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type WW(WWSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type VV(VVSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type W3(W3SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type V3(V3SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type W4(W4SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type V4(V4SEXP);
    Rcpp::traits::input_parameter< double >::type ll(llSEXP);
    Rcpp::traits::input_parameter< double >::type sum_log_factorial_Y(sum_log_factorial_YSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gradient(gradientSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type parameters(parametersSEXP);
    Rcpp::traits::input_parameter< int >::type ind(indSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type down(downSEXP);
    rcpp_result_gen = Rcpp::wrap(Mix_select_stepsize_for_a_parameter(Y, W, V, WW, VV, W3, V3, W4, V4, ll, sum_log_factorial_Y, gradient, parameters, ind, gamma, n, down));
    return rcpp_result_gen;
END_RCPP
}
// Mix_gradient_descent_for_individual_sample
Rcpp::List Mix_gradient_descent_for_individual_sample(arma::vec Y, arma::vec W, arma::vec V, double a0, double a1, double a2, double a3, double a4, double b1, double b2, double b3, double b4, double psi, double gamma, int steps, double down);
RcppExport SEXP _Imputation_Mix_gradient_descent_for_individual_sample(SEXP YSEXP, SEXP WSEXP, SEXP VSEXP, SEXP a0SEXP, SEXP a1SEXP, SEXP a2SEXP, SEXP a3SEXP, SEXP a4SEXP, SEXP b1SEXP, SEXP b2SEXP, SEXP b3SEXP, SEXP b4SEXP, SEXP psiSEXP, SEXP gammaSEXP, SEXP stepsSEXP, SEXP downSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type V(VSEXP);
    Rcpp::traits::input_parameter< double >::type a0(a0SEXP);
    Rcpp::traits::input_parameter< double >::type a1(a1SEXP);
    Rcpp::traits::input_parameter< double >::type a2(a2SEXP);
    Rcpp::traits::input_parameter< double >::type a3(a3SEXP);
    Rcpp::traits::input_parameter< double >::type a4(a4SEXP);
    Rcpp::traits::input_parameter< double >::type b1(b1SEXP);
    Rcpp::traits::input_parameter< double >::type b2(b2SEXP);
    Rcpp::traits::input_parameter< double >::type b3(b3SEXP);
    Rcpp::traits::input_parameter< double >::type b4(b4SEXP);
    Rcpp::traits::input_parameter< double >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< int >::type steps(stepsSEXP);
    Rcpp::traits::input_parameter< double >::type down(downSEXP);
    rcpp_result_gen = Rcpp::wrap(Mix_gradient_descent_for_individual_sample(Y, W, V, a0, a1, a2, a3, a4, b1, b2, b3, b4, psi, gamma, steps, down));
    return rcpp_result_gen;
END_RCPP
}
// Predict_for_individual_sample
arma::vec Predict_for_individual_sample(arma::vec W, arma::vec V, double a0, double a1, double a2, double a3, double a4, double b1, double b2, double b3, double b4);
RcppExport SEXP _Imputation_Predict_for_individual_sample(SEXP WSEXP, SEXP VSEXP, SEXP a0SEXP, SEXP a1SEXP, SEXP a2SEXP, SEXP a3SEXP, SEXP a4SEXP, SEXP b1SEXP, SEXP b2SEXP, SEXP b3SEXP, SEXP b4SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type V(VSEXP);
    Rcpp::traits::input_parameter< double >::type a0(a0SEXP);
    Rcpp::traits::input_parameter< double >::type a1(a1SEXP);
    Rcpp::traits::input_parameter< double >::type a2(a2SEXP);
    Rcpp::traits::input_parameter< double >::type a3(a3SEXP);
    Rcpp::traits::input_parameter< double >::type a4(a4SEXP);
    Rcpp::traits::input_parameter< double >::type b1(b1SEXP);
    Rcpp::traits::input_parameter< double >::type b2(b2SEXP);
    Rcpp::traits::input_parameter< double >::type b3(b3SEXP);
    Rcpp::traits::input_parameter< double >::type b4(b4SEXP);
    rcpp_result_gen = Rcpp::wrap(Predict_for_individual_sample(W, V, a0, a1, a2, a3, a4, b1, b2, b3, b4));
    return rcpp_result_gen;
END_RCPP
}
// reweighting_sum_C
arma::vec reweighting_sum_C(arma::mat Ymat, arma::mat Yflagmat, arma::vec Y, arma::vec Yflag, arma::vec prior_weight, bool ImputeAll);
RcppExport SEXP _Imputation_reweighting_sum_C(SEXP YmatSEXP, SEXP YflagmatSEXP, SEXP YSEXP, SEXP YflagSEXP, SEXP prior_weightSEXP, SEXP ImputeAllSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Ymat(YmatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Yflagmat(YflagmatSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Yflag(YflagSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type prior_weight(prior_weightSEXP);
    Rcpp::traits::input_parameter< bool >::type ImputeAll(ImputeAllSEXP);
    rcpp_result_gen = Rcpp::wrap(reweighting_sum_C(Ymat, Yflagmat, Y, Yflag, prior_weight, ImputeAll));
    return rcpp_result_gen;
END_RCPP
}
// reweighting_C
arma::vec reweighting_C(arma::mat Ymat, arma::mat Yflagmat, arma::vec Y, arma::vec Yflag);
RcppExport SEXP _Imputation_reweighting_C(SEXP YmatSEXP, SEXP YflagmatSEXP, SEXP YSEXP, SEXP YflagSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Ymat(YmatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Yflagmat(YflagmatSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Yflag(YflagSEXP);
    rcpp_result_gen = Rcpp::wrap(reweighting_C(Ymat, Yflagmat, Y, Yflag));
    return rcpp_result_gen;
END_RCPP
}
// imputation_by_samples
Rcpp::List imputation_by_samples(arma::mat data, arma::mat selected_logxx, arma::mat logxx, arma::mat zero_matrix, int n, int p, bool minbool);
RcppExport SEXP _Imputation_imputation_by_samples(SEXP dataSEXP, SEXP selected_logxxSEXP, SEXP logxxSEXP, SEXP zero_matrixSEXP, SEXP nSEXP, SEXP pSEXP, SEXP minboolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type selected_logxx(selected_logxxSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type logxx(logxxSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type zero_matrix(zero_matrixSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< bool >::type minbool(minboolSEXP);
    rcpp_result_gen = Rcpp::wrap(imputation_by_samples(data, selected_logxx, logxx, zero_matrix, n, p, minbool));
    return rcpp_result_gen;
END_RCPP
}
// imputation_by_genes
Rcpp::List imputation_by_genes(arma::mat imputed, arma::mat data_copy2, arma::uvec which_flag);
RcppExport SEXP _Imputation_imputation_by_genes(SEXP imputedSEXP, SEXP data_copy2SEXP, SEXP which_flagSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type imputed(imputedSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type data_copy2(data_copy2SEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type which_flag(which_flagSEXP);
    rcpp_result_gen = Rcpp::wrap(imputation_by_genes(imputed, data_copy2, which_flag));
    return rcpp_result_gen;
END_RCPP
}
// reweighting_with_bulk_C
arma::vec reweighting_with_bulk_C(arma::mat Ymat, arma::mat Yflagmat, arma::vec meanY, arma::vec sdY, arma::vec prior_weight);
RcppExport SEXP _Imputation_reweighting_with_bulk_C(SEXP YmatSEXP, SEXP YflagmatSEXP, SEXP meanYSEXP, SEXP sdYSEXP, SEXP prior_weightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Ymat(YmatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Yflagmat(YflagmatSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type meanY(meanYSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sdY(sdYSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type prior_weight(prior_weightSEXP);
    rcpp_result_gen = Rcpp::wrap(reweighting_with_bulk_C(Ymat, Yflagmat, meanY, sdY, prior_weight));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Imputation_calculate_weights", (DL_FUNC) &_Imputation_calculate_weights, 2},
    {"_Imputation_fitting_lasso", (DL_FUNC) &_Imputation_fitting_lasso, 3},
    {"_Imputation_log_factorial", (DL_FUNC) &_Imputation_log_factorial, 1},
    {"_Imputation_log_factorial_calculated", (DL_FUNC) &_Imputation_log_factorial_calculated, 1},
    {"_Imputation_Mix_gradient_and_LogLikelihood_for_individual_sample", (DL_FUNC) &_Imputation_Mix_gradient_and_LogLikelihood_for_individual_sample, 29},
    {"_Imputation_Mix_LogLikelihood_for_individual_sample", (DL_FUNC) &_Imputation_Mix_LogLikelihood_for_individual_sample, 21},
    {"_Imputation_Mix_select_stepsize_for_a_parameter", (DL_FUNC) &_Imputation_Mix_select_stepsize_for_a_parameter, 17},
    {"_Imputation_Mix_gradient_descent_for_individual_sample", (DL_FUNC) &_Imputation_Mix_gradient_descent_for_individual_sample, 16},
    {"_Imputation_Predict_for_individual_sample", (DL_FUNC) &_Imputation_Predict_for_individual_sample, 11},
    {"_Imputation_reweighting_sum_C", (DL_FUNC) &_Imputation_reweighting_sum_C, 6},
    {"_Imputation_reweighting_C", (DL_FUNC) &_Imputation_reweighting_C, 4},
    {"_Imputation_imputation_by_samples", (DL_FUNC) &_Imputation_imputation_by_samples, 7},
    {"_Imputation_imputation_by_genes", (DL_FUNC) &_Imputation_imputation_by_genes, 3},
    {"_Imputation_reweighting_with_bulk_C", (DL_FUNC) &_Imputation_reweighting_with_bulk_C, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_Imputation(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
