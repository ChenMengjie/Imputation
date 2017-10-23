#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
using namespace std;


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::vec calculate_weights(arma::vec z, arma::mat X){
  Environment myEnv("package:Imputation");
  Function calculate_weights_fun = myEnv["calculate_weights"];
  Rcpp::NumericVector calculate_weights_res = wrap(calculate_weights_fun(z, X));
  return calculate_weights_res;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List fitting_lasso(arma::vec y, arma::mat X, bool min){
  Environment myEnv("package:Imputation");
  Function fitting_lasso_fun = myEnv["fitting_lasso"];
  if(min){
    Rcpp::List fitting_lasso_res = wrap(fitting_lasso_fun(y, X, "min"));
    return fitting_lasso_res;
  } else {
    Rcpp::List fitting_lasso_res = wrap(fitting_lasso_fun(y, X, "1se"));
    return fitting_lasso_res;
  }
}

// [[Rcpp::export]]
double log_factorial(int Y){
  double res = 0;
  for(int kk = 1; kk <= Y; ++kk){
    res += log(kk);
  }
  return res;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::vec log_factorial_calculated(int N){

  arma::vec values = arma::zeros<arma::vec>(N+1);

  for(int kk = 1; kk <= N; ++kk){
    values(kk) = values(kk-1) + log(kk);
  }

  return values;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List Mix_gradient_and_LogLikelihood_for_individual_sample(arma::vec Y, arma::vec W, arma::vec V, arma::vec WY, arma::vec WWY, arma::vec W3Y, arma::vec W4Y,
                                                                arma::vec VY, arma::vec VVY, arma::vec V3Y, arma::vec V4Y,
                                                                arma::vec WW, arma::vec W3, arma::vec W4, arma::vec VV, arma::vec V3, arma::vec V4,
                                                                double a0, double a1, double a2, double a3, double a4,
                                                                double b1, double b2, double b3, double b4, double psi, int n, double sum_log_factorial_Y){

  arma::vec gradient = arma::zeros<arma::vec>(10);
  double Likelihood = 0;
  double common_term = -lgamma(psi) + psi*log(psi);
  double common_term2 = log(psi) - R::digamma(psi);

  for(int i = 0; i < n; ++i){

    double A_i = a0 + a1*W(i) + a2*WW(i) + a3*W3(i) + a4*W4(i);
    double B_i = b1*V(i) + b2*VV(i) + b3*V3(i) + b4*V4(i);
    double exp_Ai = exp(A_i);
    double pr =  exp_Ai*B_i;
    double ll =  pr + psi;
    double psi_Yi = psi + Y(i);
    Likelihood += lgamma(psi_Yi) - psi_Yi*log(ll) + Y(i)*(A_i + log(B_i));

    double term1 = -pr*psi_Yi/ll;
    gradient(0) += term1 + Y(i);
    gradient(1) += term1*W(i) + WY(i);
    gradient(2) += term1*WW(i) + WWY(i);
    gradient(3) += term1*W3(i) + W3Y(i);
    gradient(4) += term1*W4(i) + W4Y(i);

    double term2 = 1/B_i;
    gradient(5) += (term1*V(i) + VY(i))*term2;
    gradient(6) += (term1*VV(i) + VVY(i))*term2;
    gradient(7) += (term1*V3(i) + V3Y(i))*term2;
    gradient(8) += (term1*V4(i) + V4Y(i))*term2;
    gradient(9) += - psi_Yi/ll - log(ll) + R::digamma(psi_Yi);
  }

  gradient(9) += n*common_term2 + n;
  Likelihood -= sum_log_factorial_Y;
  Likelihood += n*common_term;
  return Rcpp::List::create(Rcpp::Named("gradient") = gradient,
                            Rcpp::Named("Likelihood") = Likelihood);

}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]


double Mix_LogLikelihood_for_individual_sample(arma::vec Y, arma::vec W, arma::vec V, arma::vec WW, arma::vec VV,
                                               arma::vec W3, arma::vec V3, arma::vec W4, arma::vec V4,
                                               double a0, double a1, double a2, double a3, double a4,
                                               double b1, double b2, double b3, double b4, double psi,
                                               int n, double sum_log_factorial_Y){

  double Likelihood = 0;
  double common_term = -lgamma(psi) + psi*log(psi);
  for(int i = 0; i < n; ++i){

    double A_i = a0 + a1*W(i) + a2*WW(i) + a3*W3(i) + a4*W4(i);
    double B_i = b1*V(i) + b2*VV(i) + b3*V3(i) + b4*V4(i);
    double exp_Ai = exp(A_i);
    double pr =  exp_Ai*B_i;
    double ll =  pr + psi;
    double psi_Yi = psi + Y(i);
    Likelihood += lgamma(psi_Yi) - psi_Yi*log(ll) + Y(i)*(A_i + log(B_i));

  }

  Likelihood -= sum_log_factorial_Y;
  Likelihood += n*common_term;
  return Likelihood;
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double Mix_select_stepsize_for_a_parameter(arma::vec Y, arma::vec W, arma::vec V, arma::vec WW, arma::vec VV,
                                           arma::vec W3, arma::vec V3, arma::vec W4, arma::vec V4, double ll, double sum_log_factorial_Y,
                                           arma::vec gradient, arma::vec parameters, int ind, double gamma, int n, double down){

  //if(ind < 5) gamma = 0.75;
  double gra = gradient(ind);
  double gra_2 = gra*gra*gamma;
  double para = parameters(ind);
  double start = sqrt(abs(para/gra))/5;
  //if(ind < 5) start = start/10;
  double a0 = parameters(0);
  double a1 = parameters(1);
  double a2 = parameters(2);
  double a3 = parameters(3);
  double a4 = parameters(4);
  double b1 = parameters(5);
  double b2 = parameters(6);
  double b3 = parameters(7);
  double b4 = parameters(8);
  double psi = parameters(9);

  double aa = start;
  double selected = para;
  double ll_prime = ll;
  while(aa > 0){
    double aa2 = aa*aa;
    double para_prime = para + aa2*gra;
    if(ind == 0){
      ll_prime = Mix_LogLikelihood_for_individual_sample(Y, W, V, WW, VV, W3, V3, W4, V4, para_prime, a1, a2, a3, a4, b1, b2, b3, b4, psi, n, sum_log_factorial_Y);
    }
    if(ind == 1){
      ll_prime = Mix_LogLikelihood_for_individual_sample(Y, W, V, WW, VV, W3, V3, W4, V4, a0, para_prime, a2, a3, a4, b1, b2, b3, b4, psi, n, sum_log_factorial_Y);
    }
    if(ind == 2){
      ll_prime = Mix_LogLikelihood_for_individual_sample(Y, W, V, WW, VV, W3, V3, W4, V4, a0, a1, para_prime, a3, a4, b1, b2, b3, b4, psi, n, sum_log_factorial_Y);
    }
    if(ind == 3){
      ll_prime = Mix_LogLikelihood_for_individual_sample(Y, W, V, WW, VV, W3, V3, W4, V4, a0, a1, a2, para_prime, a4, b1, b2, b3, b4, psi, n, sum_log_factorial_Y);
    }
    if(ind == 4){
      ll_prime = Mix_LogLikelihood_for_individual_sample(Y, W, V, WW, VV, W3, V3, W4, V4, a0, a1, a2, a3, para_prime, b1, b2, b3, b4, psi, n, sum_log_factorial_Y);
    }
    if(ind == 5){
      ll_prime = Mix_LogLikelihood_for_individual_sample(Y, W, V, WW, VV, W3, V3, W4, V4, a0, a1, a2, a3, a4, para_prime, b2, b3, b4, psi, n, sum_log_factorial_Y);
    }
    if(ind == 6){
      ll_prime = Mix_LogLikelihood_for_individual_sample(Y, W, V, WW, VV, W3, V3, W4, V4, a0, a1, a2, a3, a4, b1, para_prime, b3, b4, psi, n, sum_log_factorial_Y);
    }
    if(ind == 7){
      ll_prime = Mix_LogLikelihood_for_individual_sample(Y, W, V, WW, VV, W3, V3, W4, V4, a0, a1, a2, a3, a4, b1, b2, para_prime, b4, psi, n, sum_log_factorial_Y);
    }
    if(ind == 8){
      ll_prime = Mix_LogLikelihood_for_individual_sample(Y, W, V, WW, VV, W3, V3, W4, V4, a0, a1, a2, a3, a4, b1, b2, b3, para_prime, psi, n, sum_log_factorial_Y);
    }
    if(ind == 9){
      ll_prime = Mix_LogLikelihood_for_individual_sample(Y, W, V, WW, VV, W3, V3, W4, V4, a0, a1, a2, a3, a4, b1, b2, b3, b4, para_prime, n, sum_log_factorial_Y);
    }
    if(ll_prime - ll - aa2*gra_2 > 0 ) { //| abs(ll_prime - ll - aa2*gra_2) < 0.0001) {
      selected = para_prime;
      break;
    }
    aa = aa - start*down;
  }

  return selected;
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List Mix_gradient_descent_for_individual_sample(arma::vec Y, arma::vec W, arma::vec V, double a0, double a1, double a2, double a3, double a4,
                                                      double b1, double b2, double b3, double b4, double psi, double gamma, int steps, double down){

  int n = Y.n_elem;

  arma::vec calculated_values = log_factorial_calculated(Y.max());
  arma::vec log_factorial_Y = arma::zeros<arma::vec>(n);
  for(int i = 0; i < n; ++i){
    log_factorial_Y(i) = calculated_values(Y(i));
  }

  double sum_log_factorial_Y = sum(log_factorial_Y);

  arma::vec WW = W%W;
  arma::vec W3 = W%WW;
  arma::vec W4 = W%W3;
  arma::vec WY = W%Y;
  arma::vec WWY = W%WY;
  arma::vec W3Y = W%WWY;
  arma::vec W4Y = W%W3Y;

  arma::vec VV = V%V;
  arma::vec V3 = V%VV;
  arma::vec V4 = V%V3;
  arma::vec VY = V%Y;
  arma::vec VVY = V%VY;
  arma::vec V3Y = V%VVY;
  arma::vec V4Y = V%V3Y;
  Rcpp::List res = Mix_gradient_and_LogLikelihood_for_individual_sample(Y, W, V, WY, WWY, W3Y, W4Y, VY, VVY, V3Y, V4Y,
                                                                        WW, W3, W4, VV, V3, V4, a0, a1, a2, a3, a4, b1, b2, b3, b4, psi, n, sum_log_factorial_Y);
  arma::vec gradient = res["gradient"];
  double ll = res["Likelihood"];

  arma::vec parameters = arma::zeros<arma::vec>(10);
  parameters(0) = a0;
  parameters(1) = a1;
  parameters(2) = a2;
  parameters(3) = a3;
  parameters(4) = a4;
  parameters(5) = b1;
  parameters(6) = b2;
  parameters(7) = b3;
  parameters(8) = b4;
  parameters(9) = psi;

  double a0_prime = 0; double a1_prime = 0; double a2_prime = 0; double a3_prime = 0; double a4_prime = 0;
  double b1_prime = 0; double b2_prime = 0; double b3_prime = 0; double b4_prime = 0; double psi_prime = 1;

  for(int i = 0; i < steps; ++i){
    if(abs(gradient(0)) >= 0.0001){
      a0_prime = Mix_select_stepsize_for_a_parameter(Y, W, V, WW, VV, W3, V3, W4, V4, ll, sum_log_factorial_Y, gradient, parameters, 0, gamma, n, down);
    } else {
      a0_prime = a0;
    }
    if(abs(gradient(1)) >= 0.0001){
      a1_prime = Mix_select_stepsize_for_a_parameter(Y, W, V, WW, VV, W3, V3, W4, V4, ll, sum_log_factorial_Y, gradient, parameters, 1, gamma, n, down);
    } else {
      a1_prime = a1;
    }
    if(abs(gradient(2)) >= 0.0001){
      a2_prime = Mix_select_stepsize_for_a_parameter(Y, W, V, WW, VV, W3, V3, W4, V4, ll, sum_log_factorial_Y,  gradient, parameters, 2, gamma, n, down);
    } else {
      a2_prime = a2;
    }
    if(abs(gradient(3)) >= 0.0001){
      a3_prime = Mix_select_stepsize_for_a_parameter(Y, W, V, WW, VV, W3, V3, W4, V4, ll, sum_log_factorial_Y, gradient, parameters, 3, gamma, n, down);
    } else {
      a3_prime = a3;
    }
    if(abs(gradient(4)) >= 0.0001){
      a4_prime = Mix_select_stepsize_for_a_parameter(Y, W, V, WW, VV, W3, V3, W4, V4, ll, sum_log_factorial_Y, gradient, parameters, 4, gamma, n, down);
    } else {
      a4_prime = a4;
    }
    if(abs(gradient(5)) >= 0.0001){
      b1_prime = Mix_select_stepsize_for_a_parameter(Y, W, V, WW, VV, W3, V3, W4, V4, ll, sum_log_factorial_Y, gradient, parameters, 5, gamma, n, down);
    } else {
      b1_prime = b1;
    }
    if(abs(gradient(6)) >= 0.0001){
      b2_prime = Mix_select_stepsize_for_a_parameter(Y, W, V, WW, VV, W3, V3, W4, V4, ll, sum_log_factorial_Y,  gradient, parameters, 6, gamma, n, down);
    } else {
      b2_prime = b2;
    }
    if(abs(gradient(7)) >= 0.0001){
      b3_prime = Mix_select_stepsize_for_a_parameter(Y, W, V, WW, VV, W3, V3, W4, V4, ll, sum_log_factorial_Y, gradient, parameters, 7, gamma, n, down);
    } else {
      b3_prime = b3;
    }
    if(abs(gradient(8)) >= 0.0001){
      b4_prime = Mix_select_stepsize_for_a_parameter(Y, W, V, WW, VV, W3, V3, W4, V4, ll, sum_log_factorial_Y, gradient, parameters, 8, gamma, n, down);
    } else {
      b4_prime = b4;
    }
    if(abs(gradient(9)) >= 0.0001){
      psi_prime = Mix_select_stepsize_for_a_parameter(Y, W, V, WW, VV, W3, V3, W4, V4, ll, sum_log_factorial_Y, gradient, parameters, 9, gamma, n, down);
    } else {
      psi_prime = psi;
    }

    a0 = a0_prime; a1 = a1_prime; a2 = a2_prime; a3 = a3_prime; a4 = a4_prime;
    b1 = b1_prime; b2 = b2_prime; b3 = b3_prime; b4 = b4_prime; psi = psi_prime;

    Rcpp::List res = Mix_gradient_and_LogLikelihood_for_individual_sample(Y, W, V, WY, WWY, W3Y, W4Y, VY, VVY, V3Y, V4Y,
                                                                          WW, W3, W4, VV, V3, V4, a0, a1, a2, a3, a4, b1, b2, b3, b4, psi, n, sum_log_factorial_Y);
    arma::vec gradient_1 = res["gradient"];
    gradient = gradient_1;
    ll = res["Likelihood"];

    parameters(0) = a0;
    parameters(1) = a1;
    parameters(2) = a2;
    parameters(3) = a3;
    parameters(4) = a4;
    parameters(5) = b1;
    parameters(6) = b2;
    parameters(7) = b3;
    parameters(8) = b4;
    parameters(9) = psi;

    // Rcpp::Rcout << gradient << std::endl;
  }

  arma::vec corrected = Y;
  for(int i = 0; i < n; ++i){

    double A_i = a0 + a1*W(i) + a2*WW(i) + a3*W3(i) + a4*W4(i);
    double B_i = b1*V(i) + b2*VV(i) + b3*V3(i) + b4*V4(i);
    corrected(i) = B_i*exp(A_i);
  }

  return Rcpp::List::create(Rcpp::Named("parameters") = parameters,
                            Rcpp::Named("corrected") = corrected);
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::vec Predict_for_individual_sample(arma::vec W, arma::vec V, double a0, double a1, double a2, double a3, double a4,
                                                      double b1, double b2, double b3, double b4){

  arma::vec WW = W%W;
  arma::vec W3 = W%WW;
  arma::vec W4 = W%W3;
  arma::vec VV = V%V;
  arma::vec V3 = V%VV;
  arma::vec V4 = V%V3;

  arma::vec corrected = W;
  int n = W.n_elem;
  for(int i = 0; i < n; ++i){

    double A_i = a0 + a1*W(i) + a2*WW(i) + a3*W3(i) + a4*W4(i);
    double B_i = b1*V(i) + b2*VV(i) + b3*V3(i) + b4*V4(i);
    corrected(i) = B_i*exp(A_i);
  }

  return(corrected);

}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::vec reweighting_sum_C(arma::mat Ymat, arma::mat Yflagmat, arma::vec Y, arma::vec Yflag, arma::vec prior_weight, bool ImputeAll){

  int p = Ymat.n_cols;
  int k = Ymat.n_rows;
  double k_cut = 0.9*k;

  arma::vec res = Y;
  arma::mat Ymat_prod = Ymat%Yflagmat;

  if(ImputeAll)  {


    for(int i = 0; i < p; ++i){

      arma::vec Ymat_weights = arma::ones<arma::vec>(k);

      arma::uvec Y_ind_zero = arma::find(Yflagmat.col(i) == 0);
      double zero_num = Y_ind_zero.n_elem;
      double nonzero_num = k - zero_num;

      if(zero_num < k_cut & zero_num != 0){

        double mean_Y = sum(Ymat_prod.col(i))/nonzero_num;
        double var_Y = 1;
        if(nonzero_num != 1){
          var_Y = sum(Ymat_prod.col(i)%Ymat_prod.col(i))/nonzero_num - mean_Y*mean_Y;
          var_Y *= nonzero_num/(nonzero_num-1);
        }
        if(var_Y <= 0.01) {
          var_Y = 1;
        }
        double pN = R::dnorm(0, mean_Y, sqrt(var_Y), FALSE);
        double non_dropout = nonzero_num*pN/(zero_num + nonzero_num*pN);
        for(int j = 0; j < zero_num; ++j){
          int zero_id = Y_ind_zero(j);
          Ymat_weights(zero_id) = non_dropout;
        }
        Ymat_weights %= prior_weight;
        Ymat_weights /= sum(Ymat_weights);

      } else {

        Ymat_weights = prior_weight;

      }
      res(i) = sum(Ymat_weights%Ymat.col(i));
    }

  } else {

    arma::uvec to_impute = arma::find(Yflag == 0);
    double to_impute_num = to_impute.n_elem;
    //
    for(int i = 0; i < to_impute_num; ++i){

      arma::vec Ymat_weights = prior_weight;

      int id = to_impute(i);

      arma::uvec Y_ind_zero = arma::find(Yflagmat.col(id) == 0);
      double zero_num = Y_ind_zero.n_elem;
      double nonzero_num = k - zero_num;

      if(zero_num < k_cut & zero_num != 0){

        double mean_Y = sum(Ymat_prod.col(id))/nonzero_num;
        double var_Y = 1;
        if(nonzero_num != 1){
          var_Y = sum(Ymat_prod.col(i)%Ymat_prod.col(i))/nonzero_num - mean_Y*mean_Y;
          var_Y *= nonzero_num/(nonzero_num-1);
        }
        if(var_Y == 0) {
          var_Y = 1;
        }
        double pN = R::dnorm(0, mean_Y, sqrt(var_Y), FALSE);
        double non_dropout = nonzero_num*pN/(zero_num + nonzero_num*pN);

        for(int j = 0; j < zero_num; ++j){
          int zero_id = Y_ind_zero(j);
          Ymat_weights(zero_id) = prior_weight(zero_id)*non_dropout;
        }

        Ymat_weights = Ymat_weights/sum(Ymat_weights);
      }
      res(id) = sum(Ymat_weights%Ymat.col(id));

    }

  }
  return res;

}




// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::vec reweighting_C(arma::mat Ymat, arma::mat Yflagmat, arma::vec Y, arma::vec Yflag){

  int p = Ymat.n_cols;
  int k = Ymat.n_rows;
  double k_cut = 0.9*k;

  arma::vec res = Y;
  arma::mat Ymat_prod = Ymat%Yflagmat;

  for(int i = 0; i < p; ++i){

    arma::vec Ymat_weights = arma::ones<arma::vec>(k+1);
    arma::uvec Y_ind_zero = arma::find(Yflagmat.col(i) == 0);
    double zero_num = Y_ind_zero.n_elem;
    double nonzero_num = k - zero_num;
    if(Y(i) == 0) nonzero_num = nonzero_num + 1;

    if(zero_num < k_cut & zero_num != 0){

      double mean_Y = (sum(Ymat_prod.col(i)) + Y(i))/nonzero_num;
      double pN = R::dpois(0, mean_Y, FALSE);
      double non_dropout = nonzero_num*pN/(zero_num + nonzero_num*pN);
      for(int j = 0; j < zero_num; ++j){
        int zero_id = Y_ind_zero(j);
        Ymat_weights(zero_id) = non_dropout;
      }
      if(Y(i) == 0) {
        Ymat_weights(k) = non_dropout;
      }
    }
    Ymat_weights /= sum(Ymat_weights);
    arma::vec Values = arma::ones<arma::vec>(k+1);
    Values(arma::span(0, k-1)) = Ymat.col(i);
    Values(k) = Y(i);
    res(i) = sum(Ymat_weights%Values);
  }

  return res*(k+1);

}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List imputation_by_samples(arma::mat data, arma::mat selected_logxx, arma::mat logxx, arma::mat zero_matrix, int n, int p, bool minbool){

  arma::mat imputed = logxx;
  arma::mat t_logxx = logxx.t();
  arma::mat t_zero_matrix = zero_matrix.t();
  arma::mat sample_weights = arma::zeros<arma::mat>(n, n);

  arma::uvec ind_list = arma::zeros<arma::uvec>(n-1);
  for(int j = 0; j < n-1; ++j){
    ind_list(j) = j;
  }

  for(int j = 0; j < n; ++j){
    arma::uvec ind_new = ind_list;
    for(int l = j; l < n-1; ++l){
      ind_new(l) += 1;
    }

    Rcpp::List res = fitting_lasso(data.col(j), data.cols(ind_new), minbool);
    arma::vec coeff = res["coeff"];
    arma::uvec selected = res["selected"] ;
    selected = selected - 1;
   // Rcpp::Rcout << j << std::endl;
   //  Rcpp::Rcout << selected.n_elem << std::endl;
    if(selected.n_elem < 3){

      sample_weights.row(j).fill(-1);

    } else {

      arma::mat selected_submat = selected_logxx.cols(ind_new);
      arma::vec prior_weight = calculate_weights(selected_logxx.col(j), selected_submat.cols(selected));
      arma::uvec nonzero_ind = arma::find(prior_weight >= 0.0001);
      arma::uvec sub_selected = selected(nonzero_ind);
      arma::uvec new_lab = ind_new(sub_selected);
      arma::vec sub_prior_weight = prior_weight(nonzero_ind);
      sub_prior_weight /= sum(sub_prior_weight);
      for(int s = 0; s < new_lab.n_elem; ++s){
        sample_weights(j, new_lab(s)) = sub_prior_weight(s);
      }

      arma::mat Ymat = t_logxx.rows(new_lab);
      arma::mat Yflagmat = t_zero_matrix.rows(new_lab);

      imputed.col(j) = reweighting_sum_C(Ymat, Yflagmat, logxx.col(j), zero_matrix.col(j), sub_prior_weight, TRUE);

    }

  }

  return Rcpp::List::create(Rcpp::Named("imputed") = imputed,
                            Rcpp::Named("sample_weights") = sample_weights);
}




// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

Rcpp::List imputation_by_genes(arma::mat imputed, arma::mat data_copy2, arma::uvec which_flag){

    arma::mat imputed_gene = imputed;
    int p = data_copy2.n_cols;

    arma::vec gene_weights_1 = arma::zeros<arma::vec>(10000);
    arma::uvec gene_weights_2 = arma::zeros<arma::uvec>(10000);
    arma::vec gene_weights_3 = arma::zeros<arma::vec>(10000);

    int counter = 0; int counter_all = 0;
    arma::uvec ind_list = arma::zeros<arma::uvec>(p-1);
    for(int i = 0; i < p-1; ++i){
      ind_list(i) = i;
    }

    for(int i = 0; i < p; ++i){

      arma::uvec ind_new = ind_list;
      for(int l = i; l < p-1; ++l){
        ind_new(l) += 1;
      }
      arma::mat selected_submat = data_copy2.cols(ind_new);
      arma::mat res = arma::cor(data_copy2.col(i), selected_submat);
      arma::uvec selected = arma::find(res.row(0) > 0.8);

      int len = selected.n_elem;

      if(len == 1){
        imputed_gene.row(i) = imputed.row(selected(0));
      }

      if(len > 1){

        arma::vec prior_weight = calculate_weights(data_copy2.col(i), selected_submat.cols(selected));
        arma::uvec nonzero_ind = arma::find(prior_weight >= 0.0001);
        arma::uvec sub_selected = selected(nonzero_ind);

        if(sub_selected.n_elem > 1){

          arma::uvec new_lab = ind_new(sub_selected);
          arma::vec sub_prior_weight = prior_weight(nonzero_ind);
          sub_prior_weight /= sum(sub_prior_weight);
          arma::mat Ymat = imputed.rows(which_flag(new_lab));
          int n1 = Ymat.n_cols;
          arma::rowvec imp_neighbors = arma::zeros<arma::rowvec>(n1);
          for(int j = 0; j < n1; ++j){
            imp_neighbors(j) = sum(Ymat.col(j)%sub_prior_weight);
          }
          imputed_gene.row(i) = imp_neighbors;
          counter_all = counter + sub_selected.n_elem;
          gene_weights_1(arma::span(counter, counter_all-1)).fill(which_flag(i)+1);
          gene_weights_2(arma::span(counter, counter_all-1)) = which_flag(new_lab)+1;
          gene_weights_3(arma::span(counter, counter_all-1)) = sub_prior_weight;
          counter = counter_all;
        } else {

          arma::vec sub_prior_weight = arma::zeros<arma::vec>(len);
          double len1 = len;
          sub_prior_weight.fill(1/len1);
          arma::mat Ymat = imputed.rows(which_flag(selected));
          int n1 = Ymat.n_cols;
          arma::rowvec imp_neighbors = arma::zeros<arma::rowvec>(n1);
          for(int j = 0; j < Ymat.n_cols; ++j){
            imp_neighbors(j) = sum(Ymat.col(j)%sub_prior_weight);
          }
          imputed_gene.row(i) = imp_neighbors;
          counter_all = counter + len;
          gene_weights_1(arma::span(counter, counter_all-1)).fill(which_flag(i)+1);
          gene_weights_2(arma::span(counter, counter_all-1)) = which_flag(selected)+1;
          gene_weights_3(arma::span(counter, counter_all-1)) = sub_prior_weight;
          counter = counter_all;

        }


      }

    }

    return Rcpp::List::create(Rcpp::Named("imputed_gene") = imputed_gene,
                              Rcpp::Named("gene_weights_1") = gene_weights_1,
                              Rcpp::Named("gene_weights_2") = gene_weights_2,
                              Rcpp::Named("gene_weights_3") = gene_weights_3);
  }




// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::vec reweighting_with_bulk_C(arma::mat Ymat, arma::mat Yflagmat, arma::vec meanY, arma::vec sdY, arma::vec prior_weight){

  int p = Ymat.n_cols;
  int k = Ymat.n_rows;
  double k_cut = 0.9*k;

  arma::vec res = meanY;
  arma::mat Ymat_prod = Ymat%Yflagmat;

    for(int i = 0; i < p; ++i){

      arma::vec Ymat_weights = arma::ones<arma::vec>(k);
      arma::uvec Y_ind_zero = arma::find(Yflagmat.col(i) == 0);
      double zero_num = Y_ind_zero.n_elem;
      double nonzero_num = k - zero_num;

      if(zero_num < k_cut & zero_num != 0){

        double mean_Y_i = meanY(i);
        double sd_Y_i = sdY(i);
        double pN = R::dnorm(0, mean_Y_i, sd_Y_i, FALSE);
        double non_dropout = nonzero_num*pN/(zero_num + nonzero_num*pN);
        for(int j = 0; j < zero_num; ++j){
          int zero_id = Y_ind_zero(j);
          Ymat_weights(zero_id) = non_dropout;
        }
        Ymat_weights %= prior_weight;
        Ymat_weights /= sum(Ymat_weights);

      } else {

        Ymat_weights = prior_weight;

      }
      res(i) = sum(Ymat_weights%Ymat.col(i));
    }

  return res;

}


