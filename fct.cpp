#include <RcppArmadillo.h>
#include <math.h>
#include <Rmath.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

// MAJ le 22/10/2020

using namespace Rcpp;
// [[Rcpp::plugins("cpp11")]]


//-------------------------------------------------------------------------------------------------


// // [[Rcpp::export]]
// double rtnorm(int n, double mu, double sigma, double A, double B){
//   
//   double sample = mu+sigma*R::qnorm(R::runif(n, R::pnorm(A, mu, sigma, true, false), R::pnorm(B, mu, sigma, true, false)));
//   
//   return(sample);
// }




//-------------------------------------------------------------------------------------------------


// [[Rcpp::export]]
double log_d_Y_cpp(arma::mat Y, arma::mat M, arma::mat U_inv, arma::mat V_inv){
  
  double res = as<double>(wrap( - arma::trace( V_inv * (Y-M).t() * U_inv * (Y-M))/2  ));
  
  return(res);
  
}

// [[Rcpp::export]]
double log_d_Y_cpp_2(arma::mat Y, arma::mat M, arma::mat U_inv, arma::mat V_inv){
  arma::mat tmp = V_inv * (Y-M).t(), YM = Y-M;
  double res = 0;
  
  for(int i=0; i<Y.n_cols; i++){
    res += as<double>(wrap( tmp.row(i) * YM.col(i) ));
  }
  // double res = as<double>(wrap( - arma::trace( V_inv * (Y-M).t() * U_inv * (Y-M))/2  ));
  
  return(-res/2);
  
}


//-------------------------------------------------------------------------------------------------

// [[Rcpp::export]]
double y_Gamma_inv_y_cpp(arma::mat y_tilde, arma::sp_mat Gamma_inv){
  
  int m = y_tilde.n_rows;
  // int c = y_tilde.n_cols;
  double res=0;
  for(int i=0; i<m;i++){
    res+= as<double>(wrap(  y_tilde.row(i) * Gamma_inv * y_tilde.row(i).t() ));
  }
  return(res);
}

// [[Rcpp::export]]
arma::sp_mat update_Gamma_inv(arma::sp_mat Gamma_inv, double rho){
  int T = Gamma_inv.n_rows;
  Gamma_inv = arma::eye(T, T);
  arma::mat diag_0_Gamma =  arma::ones(T, 1) * (1+rho*rho); 
  diag_0_Gamma(0) = diag_0_Gamma(T-1) = 1;  //arma::repmat( 1+rho*rho, T, 1);
  arma::mat diag_1_Gamma = -rho* arma::ones(T-1, 1); //arma::repmat( -rho, T-1, 1);
  Gamma_inv = (arma::diagmat(diag_1_Gamma, -1) + arma::diagmat(diag_0_Gamma, 0) + arma::diagmat(diag_1_Gamma, 1))/(1-rho*rho);
  
  return(Gamma_inv );
}

//-------------------------------------------------------------------------------------------------


// [[Rcpp::export]]
arma::vec svd_Sigma_u_cpp(arma::mat Gamma_inv, arma::mat B, arma::vec sdu, arma::mat X, double sigma2, arma::vec y_tilde){
  double c = Gamma_inv.n_rows, q = B.n_rows;
  arma::mat U_l(c, c), V_l(c, c), U_r(q, q), V_r(q, q), U(c*q, c*q);
  arma::vec d_l(c), d_r(q), d(c*q);
  arma::mat S = arma::diagmat(sdu), D_r(q, q), D_l(c, c), D(c*q, c*q);
  
  arma::svd_econ(U_l, d_l, V_l, Gamma_inv, "left");
  arma::svd_econ(U_r, d_r, V_r, B.t() * S * X.t() * X * S * B , "left");
  // arma::svd(U_l, d_l, V_l, Gamma_inv);
  // arma::svd(U_r, d_r, V_r, B.t() * S * X.t() * X * S * B);
  
  D_l = arma::diagmat(d_l);
  D_r = arma::diagmat(d_r);
    
  U = arma::kron(U_l, U_r);
  D = arma::kron(D_l, D_r)/sigma2 + arma::eye(c*q, c*q);
  arma::mat D_inv = arma::diagmat(1/ D.diag());
  
  arma::mat Chol_inv = U * arma::diagmat(sqrt(D_inv));
  
  // Rcout << " 1 \n";
  arma::vec mean_u_t_t = Chol_inv.t() * arma::kron(Gamma_inv, B.t() * S * X.t()) * y_tilde / sigma2; 
  
  
  // Rcout << Chol_inv * (arma::randn(c*q) + mean_u_t_t) <<"  \n";
  
  arma::vec u_t = kron(arma::eye(c, c), B)  * (Chol_inv * (arma::randn(c*q) + mean_u_t_t) );
  
  // List toto = List::create(Named("chol_Sigma") = Chol_inv, Named("mean") = Chol_inv*mean_u_t_t);
  // return(kron(arma::eye(c, c), B)  * Chol_inv*mean_u_t_t);
  return(u_t);
}


//-------------------------------------------------------------------------------------------------

// [[Rcpp::export]]
List svd_Sigma_u_cpp_2(arma::sp_mat Gamma_inv, arma::mat B, arma::vec sdu, arma::mat X, double sigma2, arma::vec y_tilde){
  double c = Gamma_inv.n_rows, q = B.n_rows;
  arma::mat U_l(c, c), V_l(c, c), U_r(q, q), V_r(q, q), U(c*q, c*q);
  arma::vec d_l(c), d_r(q), d(c*q);
  arma::mat S = arma::diagmat(sdu), D_r(q, q), D_l(c, c), D(c*q, c*q);
  
  // arma::svd_econ(U_l, d_l, V_l, Gamma_inv, "left");
  arma::svd_econ(U_r, d_r, V_r, B.t() * S * X.t() * X * S * B , "left");
  arma::svds(U_l, d_l, V_l, Gamma_inv, c);
  // arma::svd(U_r, d_r, V_r, B.t() * S * X.t() * X * S * B);
  
  D_l = arma::diagmat(d_l);
  D_r = arma::diagmat(d_r);
  
  U = arma::kron(U_l, U_r);
  D = arma::kron(D_l, D_r)/sigma2 + arma::eye(c*q, c*q);
  arma::mat D_inv = arma::diagmat(1/ D.diag());
  
  arma::mat Chol_inv = U * arma::diagmat(sqrt(D_inv));
  
  arma::vec mean_u_t_t = Chol_inv.t() * arma::kron(arma::conv_to<arma::mat>::from(Gamma_inv), B.t() * S * X.t()) * y_tilde / sigma2; 
  
  
  // Rcout << Chol_inv * (arma::randn(c*q) + mean_u_t_t) <<"  \n";
  
  arma::vec u_tt =  (Chol_inv * (arma::randn(c*q) + mean_u_t_t) );
  
  arma::vec u_t = kron(arma::eye(c, c), B)  * u_tt;
  List ret = List::create(Named("u_tt") = u_tt, Named("u_t") = u_t);

  return(ret);
}



//-------------------------------------------------------------------------------------------------
// [[Rcpp::export]]
arma::sp_mat kron_cpp(arma::mat A, arma::mat B){
  return(arma::sp_mat(kron(A, B)));
}


// _____________________________________________________________________________


// [[Rcpp::export]]
arma::sp_mat kron_cpp_s_s(arma::sp_mat A, arma::sp_mat B){
  return(arma::sp_mat(kron(A, B)));
}


// _____________________________________________________________________________



// [[Rcpp::export]]
arma::mat my_chol2inv(arma::mat B){
  arma::mat B_inv = inv(B);
  return(B_inv.t() * B_inv);
  // arma::mat iR_star = inv_sympd(B * B.t()) ;
  // return(iR_star);
}



// _____________________________________________________________________________




// [[Rcpp::export]]
double l_p_Y_theta_star_cpp(int niv, int q_rand, arma::mat U, arma::vec theta_vec, arma::mat iR_star){
  // iR_star <- chol2inv(t(B_star))
  double lol = as<double>(wrap(    trace(U.t() * iR_star * U)   )) ;
  double l_p_Y_theta_star = as<double>(wrap(   - (log(2*M_PI)) * niv * q_rand/2 - niv * sum(log(sin(theta_vec))) - lol /2    ));
  return(l_p_Y_theta_star);
}


// _____________________________________________________________________________



// [[Rcpp::export]]
double l_p_Y_theta_star_cpp_2(int niv, int q_rand, arma::mat U, arma::vec theta_vec, arma::mat B, int k){
  

  arma::mat iR_star = my_chol2inv(B) ;
  double lol = as<double>(wrap(    trace(U.t() * iR_star * U)   )) ;
  double l_p_Y_theta_star = as<double>(wrap(   - niv * q_rand/2*(log(2*M_PI)) +  sum(log(sin(theta_vec))) * (k - niv) - lol /2    ));
  return(l_p_Y_theta_star);
}


// [[Rcpp::export]]
List l_p_Y_theta_star_cpp_3(int niv, int q_rand, arma::mat U, arma::vec theta_vec, arma::mat B, int k){
  
  
  arma::mat iR_star = my_chol2inv(B) ;
  double lol = as<double>(wrap(    trace(U.t() * iR_star * U)   )) ;
  double l_p_Y_theta_star = as<double>(wrap(   - niv * q_rand/2*(log(2*M_PI)) +  sum(log(sin(theta_vec))) * (k - niv) - lol /2    ));
  List res = List::create(Named("log_density") = l_p_Y_theta_star, Named("iR") = iR_star);
  return(res);
}

// _____________________________________________________________________________



// // [[Rcpp::export]]
// double l_p_Y_theta_star_cpp_3(int niv, int q_rand, arma::mat U, arma::vec theta_vec, arma::mat B){
//   arma::mat iR_star = inv_sympd(B * B.t()) ;
//   double lol = as<double>(wrap(    trace(U.t() * iR_star * U)   )) ;
//   double l_p_Y_theta_star = as<double>(wrap(   - niv * q_rand/2*(log(2*M_PI)) +  sum(log(sin(theta_vec))) * (- niv) - lol /2    ));
//   return(l_p_Y_theta_star);
// }


// _____________________________________________________________________________

// [[Rcpp::export]]
arma::vec my_mu_part(arma::mat iS_part, arma::mat tmp, arma::vec y_tild){
  arma::vec mu_part = iS_part * tmp.t() * y_tild ;
  return(mu_part);
}



//  from shrinkTVP package Bitto & Frühwirth-Schnatter 2018
/*---------------------------------------------------------------------------*/
/*                                                                           */
/* Original implementation in R code (R package "Bessel" v. 0.5-3) by        */
/*   Martin Maechler, Date: 23 Nov 2009, 13:39                               */
/*                                                                           */
/* Translated into C code by Kemal Dingic, Oct. 2011.                        */
/*                                                                           */
/* Modified by Josef Leydold on Tue Nov  1 13:22:09 CET 2011                 */
/*                                                                           */
/* Translated into C++ code by Peter Knaus, Mar. 2019.                       */
/*                                                                           */
/*---------------------------------------------------------------------------*/

double unur_bessel_k_nuasympt(double x, double nu, bool islog, bool expon_scaled){
  double M_LNPI = 1.14472988584940017414342735135;

  double z;
  double sz, t, t2, eta;
  double d, u1t,u2t,u3t,u4t;
  double res;


  z = x / nu;

  sz = hypot(1,z);
  t = 1. / sz;
  t2 = t*t;

  if (expon_scaled){
    eta = (1./(z + sz));
  } else {
    eta = sz;
  }

  eta += log(z) - log1p(sz);
  u1t = (t * (3. - 5.*t2))/24.;
  u2t = t2 * (81. + t2*(-462. + t2 * 385.))/1152.;
  u3t = t*t2 * (30375. + t2 * (-369603. + t2 * (765765. - t2 * 425425.)))/414720.;
  u4t = t2*t2 * (4465125.
                   + t2 * (-94121676.
                   + t2 * (349922430.
                   + t2 * (-446185740.
                   + t2 * 185910725.)))) / 39813120.;
                   d = (-u1t + (u2t + (-u3t + u4t/nu)/nu)/nu)/nu;

                   res = log(1.+d) - nu*eta - 0.5*(log(2.*nu*sz) - M_LNPI);

                   if (islog){
                     return res;
                   } else {
                     return exp(res);
                   }
}



// _____________________________________________________________________________




//  from shrinkTVP package Bitto & Frühwirth-Schnatter 2018
double log_ratio_value_marginalBFS(int d, double proposal, double old_val, double scale_par, arma::vec param_vec, double b1, double b2) {
  arma::vec besselKvalue_partA(d, arma::fill::none);
  arma::vec besselKvalue_partB(d, arma::fill::none);
  long double par1A = std::abs(proposal - 0.5);
  long double par1B = std::abs(old_val - 0.5);
  for (int j = 0; j < d; j++){
    double par2A = std::exp(0.5*std::log(proposal) + 0.5*std::log(scale_par) + std::log(std::abs(param_vec(j))));
    double par2B = std::exp(0.5*std::log(old_val) + 0.5*std::log(scale_par) + std::log(std::abs(param_vec(j))));


    if(par1A < 50 and par2A < 50){
      besselKvalue_partA(j) = std::log(R::bessel_k(par2A, par1A, true)) - par2A;
    }else{
      besselKvalue_partA(j) = unur_bessel_k_nuasympt(par2A, par1A, true, false);
    }

    if(par1B < 50 and par2B < 50){
      besselKvalue_partB(j) = std::log(R::bessel_k(par2B, par1B, true)) - par2B;
    }else{
      besselKvalue_partB(j) = unur_bessel_k_nuasympt(par2B, par1B, true, false);
    }

  }

  //gamma prior
  double partA = (b1 - 1 + 1 + d/4.0)*(std::log(proposal) - std::log(old_val));
  double partB = (d/2.0*std::log(scale_par) - d*std::log(2) -
                  b2 + arma::as_scalar(arma::sum(arma::log(arma::abs(param_vec)))))*(proposal - old_val);
  double partC = d/2.0*(std::log(proposal)*proposal - std::log(old_val)*old_val);
  double partD =  - d*(std::lgamma(proposal + 1) - log(proposal) - std::lgamma(old_val +1) + log(old_val));
  double partE = arma::as_scalar(arma::sum(besselKvalue_partA -besselKvalue_partB));
  double res = partA + partB + partC +partD + partE;

  return res;
}

// _____________________________________________________________________________

//  from shrinkTVP package Bitto & Frühwirth-Schnatter 2018
void res_protector(double& x){
  if (std::abs(x) < DBL_MIN * std::pow(10, 10)){
    double sign = std::copysign(1, x);
    x = DBL_MIN * std::pow(10, 10) * sign;
  }
}

// _____________________________________________________________________________
//  from shrinkTVP package Bitto & Frühwirth-Schnatter 2018

// [[Rcpp::export]]
double MH_step(double current_val, double c_tuning_par, int d, double scale_par, arma::vec param_vec, double b, double nu,
               double hyp1, double hyp2){

  double b1 = nu;
  double b2 = nu * b;

  double old_value = current_val;
  double log_prop = R::rnorm(std::log(old_value), c_tuning_par);
  double proposal = std::exp(log_prop);

  double unif = R::runif(0, 1);

  double log_R = log_ratio_value_marginalBFS(d, proposal, old_value, scale_par, param_vec, b1, b2);

  double res;
  if (std::log(unif) < log_R){
    res = proposal;
  } else {
    res = old_value;
  }

  res_protector(res);

  return res;

}




