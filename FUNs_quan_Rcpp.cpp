#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector get_basis_a(NumericVector kvec, 
                          NumericVector qvec, 
                          NumericVector thetavec) {
  int L = thetavec.size();
  NumericVector out(L);
  for(int i = 0; i < L; ++i) {
    if(kvec[i] < 0.5){
      double denval = R::qnorm(kvec[i+1], 0.0, 1.0, true, false);
      out[i] = qvec[i+1] - thetavec[i]*denval;
    } else {
      double denval = R::qnorm(kvec[i], 0.0, 1.0, true, false);
      out[i] = qvec[i] - thetavec[i]*denval;
    }
  }
  return out;
}

// [[Rcpp::export]]
NumericVector get_den(NumericVector xval, 
                      NumericMatrix intervalmat, 
                      NumericVector avec, 
                      NumericVector thetavec) {
  int L = thetavec.size();
  int n = xval.size();
  NumericVector out(n);
  for(int i = 0; i < n; ++i) {
    double deni = 0;
    for(int j = 0; j < L; ++j){
      LogicalVector indexvec = LogicalVector::create(false, false);
      indexvec[0] = (intervalmat(j, 0) <= xval[i]);
      indexvec[1] = (xval[i] < intervalmat(j, 1));
      bool index = is_true(all(indexvec)); 
      double indexuse = 1;
      if(index == TRUE){
        indexuse = 1;
      }else{
        indexuse = 0;
      }
      double denval = R::dnorm(xval[i], avec[j], thetavec[j], false);
      deni += indexuse*denval;
    }
    out[i] = deni;
  }
  return out;
}

// [[Rcpp::export]]
NumericVector comp_ll_alpha(NumericVector alph_c,
                            NumericMatrix x_sim_mat,
                            NumericMatrix theta_mat_i,
                            Function f,
                            double L,
                            double shape,
                            String basis_fun,
                            double ntimes) {
  NumericVector out(ntimes);
  for(int tt = 0; tt < ntimes; ++tt){
    NumericVector xval_tt_0 = x_sim_mat(tt, _);
    NumericVector xval_tt = na_omit(xval_tt_0);
    NumericVector theta_mat_i_tt = theta_mat_i(_, tt);
    NumericVector l_curr_t_vec = f(xval_tt, L, theta_mat_i_tt, alph_c[tt], basis_fun, shape);
    NumericVector ll_curr_t_vec = log(l_curr_t_vec);
    out[tt] = sum(ll_curr_t_vec);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector get_basis_a_gamma(NumericVector kvec, 
                                NumericVector qvec, 
                                NumericVector thetavec,
                                double shape) {
  int L = thetavec.size();
  NumericVector out(L);
  for(int i = 0; i < L; ++i) {
    if(kvec[i] < 0.5){
      double denval = R::qgamma(kvec[i+1], shape, 1.0, true, false);
      out[i] = qvec[i+1] - thetavec[i]*denval;
    } else {
      double denval = R::qgamma(kvec[i], shape, 1.0, true, false);
      out[i] = qvec[i] - thetavec[i]*denval;
    }
  }
  return out;
}


// [[Rcpp::export]]
NumericVector get_den_gamma(NumericVector xval, 
                            NumericMatrix intervalmat, 
                            NumericVector avec, 
                            NumericVector thetavec,
                            double shape) {
  int L = thetavec.size();
  int n = xval.size();
  NumericVector out(n);
  for(int i = 0; i < n; ++i) {
    double deni = 0;
    for(int j = 0; j < L; ++j){
      LogicalVector indexvec = LogicalVector::create(false, false);
      indexvec[0] = (intervalmat(j, 0) <= xval[i]);
      indexvec[1] = (xval[i] < intervalmat(j, 1));
      bool index = is_true(all(indexvec)); 
      double indexuse = 1;
      if(index == TRUE){
        indexuse = 1;
      }else{
        indexuse = 0;
      }
      double denval0 = R::dgamma((xval[i]-avec[j])/thetavec[j], shape, 1, false);
      double denval = denval0/thetavec[j];
      deni += indexuse*denval;
    }
    out[i] = deni;
  }
  return out;
}


// [[Rcpp::export]]
Rcpp::List get_Mlist(arma::mat BK_beta,
                     arma::vec omega_vec,
                     List theta_pre){
  int ntotal = theta_pre.size();
  Rcpp::List out(ntotal);
  arma::mat BK_beta_t = BK_beta.t()*BK_beta;
  for(int i=0; i<ntotal; ++i){
    arma::mat theta_pre_i = theta_pre[i];
    double omega_vec_i = omega_vec[i];
    arma::mat out_i = BK_beta_t*omega_vec_i;
    arma::mat out_i1 = out_i + theta_pre_i;
    out[i] = out_i1.i();
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::List get_mulist(List theta_prod,
                      List M_list,
                      arma::vec omega_vec,
                      arma::vec z,
                      arma::mat BK_beta){
  int ntotal = theta_prod.size();
  Rcpp::List out(ntotal);
  for(int i=0; i<ntotal; ++i){
    arma::mat theta_prod_i = theta_prod[i];
    arma::mat M_i = M_list[i];
    double omega_vec_i = omega_vec[i];
    double z_i = z[i];
    arma::mat out_i = BK_beta.t()*omega_vec_i*z_i;
    arma::mat out_i1 = out_i + theta_prod_i;
    out[i] = M_i*out_i1;
  }
  return out;
}

// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

// [[Rcpp::export]]
arma::mat sample_basiscoef(List mu_list,
                           List sigma_list,
                           int ntheta){
  int ntotal = mu_list.size();
  arma::mat out;
  out.reshape(ntotal, ntheta);
  for(int i=0; i<ntotal; ++i){
    arma::vec mu_i = mu_list[i];
    arma::mat sigma_i = sigma_list[i];
    arma::mat yreturn = mvrnormArma(1, mu_i, sigma_i);
    out.row(i) = yreturn.row(0);
  }
  return out;
}


// [[Rcpp::export]]
arma::mat compq_cpp_XZonly(arma::mat beta_x_post,
                           arma::mat beta_z_post,
                           List Xdesign_quan_post,
                           arma::mat Zdesign,
                           int ntotal, int npost){
  int nbetaz = beta_z_post.n_cols;
  int nbetax = beta_x_post.n_cols;
  arma::mat out;
  out.reshape(ntotal, npost);
  for(int i=0; i<npost; ++i){
    arma::mat beta_z_i = beta_z_post.row(i);
    beta_z_i.reshape(nbetaz,1);
    arma::mat beta_x_i = beta_x_post.row(i);
    beta_x_i.reshape(nbetax,1);
    arma::mat Xdesign_i = Xdesign_quan_post[i];
    arma::mat comp1 = Zdesign*beta_z_i;
    arma::mat comp2 = Xdesign_i*beta_x_i;
    arma::mat comp = comp1 + comp2;
    arma::mat compout = 1/(1 + exp(comp));
    out.col(i) = compout;
  }
  return out;
}

// [[Rcpp::export]]
arma::mat compq_cpp(arma::mat beta_x_post,
                    arma::mat beta_z_post,
                    arma::mat rand_post,
                    List Xdesign_quan_post,
                    arma::mat Zdesign,
                    int ntotal, int npost){
  int nbetaz = beta_z_post.n_cols;
  int nbetax = beta_x_post.n_cols;
  arma::mat out;
  out.reshape(ntotal, npost);
  for(int i=0; i<npost; ++i){
    arma::mat beta_z_i = beta_z_post.row(i);
    beta_z_i.reshape(nbetaz,1);
    arma::mat beta_x_i = beta_x_post.row(i);
    beta_x_i.reshape(nbetax,1);
    arma::mat Xdesign_i = Xdesign_quan_post[i];
    arma::mat comp1 = Zdesign*beta_z_i;
    arma::mat comp2 = Xdesign_i*beta_x_i;
    arma::mat comp3 = rand_post.row(i);
    comp3.reshape(ntotal,1);
    arma::mat comp = comp1 + comp2 + comp3;
    arma::mat compout = 1/(1 + exp(comp));
    out.col(i) = compout;
  }
  return out;
}
